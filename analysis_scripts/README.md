# mad_simulations

Resample every population, then compute **MAD** — the **m**inimum number of
**a**lleles whose cumulative allele frequency first exceeds a cutoff.

Each population is randomly sampled *with replacement* so that the quantity of
interest can be studied across resampling replicates. The simulation runs
independently at two population levels — `pop1` (continental) and `pop2`
(sub-continental). The run can optionally be restricted to a subset of loci with
`--loci` (default: all loci). Three analyses can be produced:

* **MAD** (default) — reported per level, population code, locus, cutoff, and
  replicate.
* **`--allele-ci`** — the bootstrap frequency of *every* observed allele,
  reported per level, code, locus, haplotype, and replicate, with a percentile
  confidence interval summarised per allele.
* **`--carrier-mad`** — MAD defined by carrier frequency instead of allele
  frequency: greedily pick the allele carried by the most people, drop
  everyone who carries it, and repeat until the cutoff fraction of people is
  covered (a port of `carrier_mad.py`). Because coverage is about *people*,
  this mode bootstraps whole samples with replacement rather than individual
  allele observations, and needs an extra `sample` column in the input.

This is a multithreaded Rust port of the original `mad_simulations.py`. Work is
fanned out over `replicates × groups` tasks with [Rayon](https://docs.rs/rayon).
Each task seeds its own RNG deterministically from `(seed, rep, group)`, so output
is **identical regardless of thread count** — parallelism never changes results.

## Build

```sh
cargo build --release
# binary at ./target/release/mad_simulations
```

## Input format

One CSV row per observed allele, with these columns (extra columns are ignored):

| column      | meaning                                  |
|-------------|------------------------------------------|
| `haplotype` | the allele / haplotype label             |
| `pop1`      | continental population code              |
| `pop2`      | sub-continental population code          |
| `locus`     | locus identifier                         |
| `sample`    | person/sample id (required only for `--carrier-mad`) |

```csv
haplotype,pop1,pop2,locus
A,EUR,CEU,L1
B,EUR,GBR,L1
X,AFR,YRI,L2
```

`--carrier-mad` needs a `sample` column so rows carried by the same person
(e.g. both haplotypes of a diploid genotype) stay grouped together:

```csv
sample,haplotype,pop1,pop2,locus
p1,A,EUR,CEU,L1
p1,B,EUR,CEU,L1
p2,A,EUR,CEU,L1
p2,A,EUR,CEU,L1
```

Inputs may be **plain or gzip-compressed** — gzip is detected by magic bytes, so
`.csv`, `.csv.gz`, `.gz`, and even extensionless gzip streams all work, and you can
mix compressed and uncompressed files in a single invocation. Multiple input files
are concatenated.

## Usage

Exactly one resampling mode is required:

* **`--size N` / `-n N`** — resample to a fixed `N` allele observations per
  population/locus (or `N` *people* with `--carrier-mad`), so MAD is compared
  on an equal sample-size footing.
* **`--population-size` / `-N`** — bootstrap each group at its own observed allele
  count (or its own number of distinct people with `--carrier-mad`) — a classic
  bootstrap that preserves each group's native sample size.

```sh
# Fixed sample size of 50, 1000 replicates, seed 42, using 8 threads
./target/release/mad_simulations input.csv -n 50 -r 1000 -s 42 -t 8

# Population-size bootstrap over a gzipped input
./target/release/mad_simulations input.csv.gz --population-size --reps 200

# Multiple inputs (plain + gzip), custom output path
./target/release/mad_simulations a.csv b.csv.gz -n 50 -o results.csv

# Per-allele bootstrap frequency CIs instead of MAD
./target/release/mad_simulations input.csv -n 50 --allele-ci

# Carrier-frequency MAD: bootstrap 50 people per replicate (needs a `sample` column)
./target/release/mad_simulations input.csv -n 50 --carrier-mad

# Restrict to a subset of loci (comma-separated or repeated -l)
./target/release/mad_simulations input.csv -n 50 -l L1,L2
```

### Options

| flag                        | default              | description                                            |
|-----------------------------|----------------------|--------------------------------------------------------|
| `<INPUT>...`                | —                    | one or more input CSVs (plain or gzip)                 |
| `-n, --size <N>`            | —                    | fixed sample size per group (mutually exclusive with `-N`) |
| `-N, --population-size`     | —                    | bootstrap at each group's observed count               |
| `--allele-ci`               | off                  | report per-allele bootstrap frequency CIs instead of MAD (mutually exclusive with `--carrier-mad`) |
| `--carrier-mad`              | off                  | report carrier-frequency MAD instead of allele-frequency MAD; requires a `sample` column (mutually exclusive with `--allele-ci`) |
| `-l, --loci <LOCUS>...`     | all loci             | restrict the run to this subset of loci (comma-separated or repeated) |
| `-r, --reps <R>`            | `100`                | number of resampling replicates                        |
| `-s, --seed <S>`            | `0`                  | random seed (reproducible)                             |
| `-c, --ci-level <PCT>`      | `95`                 | confidence level (percent) for the reported CI bounds  |
| `-t, --threads <T>`         | all cores            | number of worker threads                               |
| `-o, --output <PATH>`       | `mad_simulations_<tag>.csv` | output CSV path                                 |

## Output

Two files are written (`<tag>` is `n<SIZE>` or `popsize`). In MAD mode the
default prefix is `mad_simulations_`; with `--allele-ci` it is
`allele_ci_simulations_`; with `--carrier-mad` it is
`carrier_mad_simulations_`. The `ci<LEVEL>` column prefix follows `--ci-level`
(e.g. `ci95` by default, `ci99` for `--ci-level 99`).

### MAD (default)

* **`mad_simulations_<tag>.csv`** — one row per replicate × group × cutoff, with
  columns: `rep, level, code, locus, cutoff, n_alleles, cumulative_pop_af,
  n_sampled_alleles, sample_size`.
* **`mad_simulations_<tag>_summary.csv`** — `n_alleles` aggregated per
  `(level, code, locus, cutoff)`: `mean, median, std, min, max, ci95_low, ci95_high`.

Cutoffs span every 0.01 from 0.50 to 0.90 inclusive (41 values).

### Allele-frequency CI (`--allele-ci`)

* **`allele_ci_simulations_<tag>.csv`** — one row per replicate × group × observed
  allele, with columns: `rep, level, code, locus, haplotype, allele_af, sample_size`.
  Alleles not drawn in a given replicate carry `allele_af = 0.0`, so each allele's
  CI reflects its full bootstrap distribution.
* **`allele_ci_simulations_<tag>_summary.csv`** — `allele_af` aggregated per
  `(level, code, locus, haplotype)`: `mean, median, std, min, max, ci95_low, ci95_high`.

### Carrier-frequency MAD (`--carrier-mad`)

* **`carrier_mad_simulations_<tag>.csv`** — one row per replicate × group ×
  cutoff, with columns: `rep, level, code, locus, cutoff, n_alleles,
  cumulative_carrier_frac, n_available_alleles, sample_size`. Here `sample_size`
  is the number of *people* bootstrapped into the replicate (not allele
  observations), and `cumulative_carrier_frac` is the fraction of them covered
  by the chosen alleles.
* **`carrier_mad_simulations_<tag>_summary.csv`** — `n_alleles` aggregated per
  `(level, code, locus, cutoff)`, same shape as the default MAD summary.

Note this mode has no `pop2` vs `pop1` interaction beyond the usual two levels
— there is no pooled `global` level (that's `--allele-ci`-only).

In all modes the summary uses sample std (`ddof=1`; blank when n < 2) and
linear-interpolated percentiles for the CI bounds. A custom `-o name.csv`
produces `name_summary.csv` alongside it.

## Reproducibility note

Results are statistically equivalent to the original Python version but **not
bit-for-bit identical** to its specific random draws: NumPy's `default_rng` (PCG64)
stream can't be reproduced by a different RNG. The procedure and summary statistics
are the same; only the exact per-replicate resamples differ. For a given seed, the
Rust version is fully reproducible against itself and across any thread count.
