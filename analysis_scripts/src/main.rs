//! Resample every population, then compute MAD.
//!
//! Each population is randomly sampled *with replacement* so that MAD can be
//! studied across resampling replicates. Two modes are supported:
//!
//! * `--size N`: resample down/up to a fixed `N` allele observations per
//!   population/locus so that MAD can be compared on an equal sample-size footing.
//! * `--population-size`: resample with replacement up to each population's own
//!   observed allele count (alleles per locus per population), i.e. a classic
//!   bootstrap that preserves the native sample size of every group.
//!
//! The simulation is run independently at each population level in `LEVELS`
//! (`pop1` continental, `pop2` sub-continental). The input is one row per
//! observed allele (columns: `haplotype`, `pop1`, `pop2`, `locus`); allele
//! frequencies are computed on the fly from the resampled alleles.
//!
//! Three analyses can be produced from the same resampling:
//!
//! * MAD (default): the minimum number of alleles whose cumulative frequency
//!   first exceeds each cutoff, reported per level, code, locus, cutoff, and
//!   replicate.
//! * `--allele-ci`: the bootstrap frequency of *every* observed allele, reported
//!   per level, code, locus, haplotype, and replicate, with a percentile
//!   confidence interval summarised per allele. In addition to `pop1`/`pop2`, an
//!   extra `global` level pools every population together and resamples each
//!   locus up to the total size of all pops: in `--population-size` mode that is
//!   the total observed allele count for the locus, and in `--size N` mode it is
//!   `N` times the number of (sub-continental) pop2 populations observed for
//!   that locus. This is computed once (not separately per population level).
//! * `--carrier-mad`: MAD defined by carrier frequency instead of allele
//!   frequency — greedily pick the allele carried by the most people, drop
//!   everyone who carries it, and repeat until the cutoff fraction of people is
//!   covered (a port of `carrier_mad.py`'s `min_alleles_for_coverage`). Because
//!   coverage is about *people*, not allele draws, this mode bootstraps whole
//!   samples with replacement (each draw brings along everything that sample
//!   carries at the locus) rather than individual allele observations, and the
//!   input needs an extra `sample` column identifying which rows belong to the
//!   same person.
//!
//! This is a multithreaded Rust port of the original `mad_simulations.py`. Work
//! is fanned out over (replicate x group) tasks with Rayon. Each task seeds its
//! own RNG deterministically from `(seed, rep, level, group)`, so results are
//! fully reproducible regardless of thread count or scheduling order. (The exact
//! random stream differs from NumPy's `default_rng`, so individual draws will not
//! match the Python version bit-for-bit, but the statistical procedure is identical.)

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::process::exit;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

/// Population levels the simulation is run over, independently.
const LEVELS: [&str; 2] = ["pop1", "pop2"];

#[derive(Parser, Debug)]
#[command(
    name = "mad_simulations",
    about = "Resample every population, then compute MAD per level, code, locus, cutoff, and replicate."
)]
struct Args {
    /// Input CSV(s) (optionally gzip-compressed), one row per allele, with
    /// columns: haplotype, pop1, pop2, locus.
    #[arg(required = true)]
    input: Vec<String>,

    /// Resample to a fixed number of allele observations per population/locus.
    #[arg(short = 'n', long, conflicts_with = "population_size")]
    size: Option<usize>,

    /// Resample with replacement up to each population's own observed allele
    /// count (alleles per locus per population).
    #[arg(short = 'N', long)]
    population_size: bool,

    /// Report bootstrap frequency confidence intervals for every observed allele
    /// (per level, code, locus, haplotype) instead of MAD.
    #[arg(long, conflicts_with = "carrier_mad")]
    allele_ci: bool,

    /// Report carrier-frequency MAD instead of allele-frequency MAD: greedily
    /// pick the allele carried by the most people, drop everyone who carries
    /// it, and repeat until the cutoff fraction of people is covered. Bootstraps
    /// whole samples with replacement rather than individual allele
    /// observations, so the input needs a `sample` column.
    #[arg(long, conflicts_with = "allele_ci")]
    carrier_mad: bool,

    /// Restrict the simulation to this subset of loci. Repeat the flag or
    /// comma-separate values (e.g. `-l L1,L2`). Default: all loci in the input.
    #[arg(short = 'l', long, value_delimiter = ',')]
    loci: Vec<String>,

    /// Number of resampling replicates.
    #[arg(short = 'r', long, default_value_t = 100)]
    reps: u64,

    /// Random seed for reproducibility.
    #[arg(short = 's', long, default_value_t = 0)]
    seed: u64,

    /// Confidence level (percent) for the reported CI, e.g. 99 for a 99% CI.
    #[arg(short = 'c', long, default_value_t = 95.0)]
    ci_level: f64,

    /// Output CSV path (default: mad_simulations_<size_tag>.csv).
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Number of worker threads (default: all available cores).
    #[arg(short = 't', long)]
    threads: Option<usize>,
}

/// Every 0.01 from 0.50 to 0.90 (inclusive): 41 cutoffs.
fn cutoffs() -> Vec<f64> {
    (0..=40).map(|i| round_half_even(0.50 + 0.01 * i as f64, 2)).collect()
}

/// Population level name for the pooled-across-all-pops group (`--allele-ci`).
const GLOBAL_LEVEL: &str = "global";

/// One (level, code, locus) group ready for resampling.
struct Group {
    level: &'static str,
    code: String,
    locus: String,
    /// One entry per observed allele, encoded as a dense allele id `0..n_alleles`.
    obs: Vec<u32>,
    /// Number of distinct alleles observed in this group.
    n_alleles: usize,
    /// Haplotype label for each allele id, indexed by id (`0..n_alleles`).
    labels: Vec<String>,
    /// For the `global` level in `--size` mode: number of pop2 populations
    /// observed for this locus (the resample size is `--size` x this). Unused
    /// (0) for the per-population levels.
    n_pops: usize,
}

/// One output row of the per-replicate simulation table.
struct Row {
    rep: u64,
    level: &'static str,
    code: String,
    locus: String,
    cutoff: f64,
    n_alleles: usize,
    cumulative_pop_af: f64,
    n_sampled_alleles: usize,
    sample_size: usize,
}

/// One output row of the per-replicate allele-frequency table (`--allele-ci`).
struct AlleleRow {
    rep: u64,
    level: &'static str,
    code: String,
    locus: String,
    haplotype: String,
    allele_af: f64,
    sample_size: usize,
}

/// One (level, code, locus) group ready for carrier-frequency resampling
/// (`--carrier-mad`). Unlike `Group`, the resampling unit here is the person,
/// not the allele observation.
struct CarrierGroup {
    level: &'static str,
    code: String,
    locus: String,
    /// One entry per distinct sample; each is the dense allele ids (`0..n_alleles`)
    /// that sample carries at this locus (e.g. 1 or 2 for a diploid genotype).
    persons: Vec<Vec<u32>>,
    /// Number of distinct alleles observed in this group.
    n_alleles: usize,
}

/// One output row of the per-replicate carrier-MAD table (`--carrier-mad`).
struct CarrierRow {
    rep: u64,
    level: &'static str,
    code: String,
    locus: String,
    cutoff: f64,
    n_alleles: usize,
    cumulative_carrier_frac: f64,
    n_available_alleles: usize,
    sample_size: usize,
}

fn main() {
    let args = Args::parse();

    if args.size.is_none() && !args.population_size {
        eprintln!("error: one of --size/-n or --population-size/-N is required");
        exit(2);
    }

    if !(args.ci_level > 0.0 && args.ci_level < 100.0) {
        eprintln!("error: --ci-level must be between 0 and 100 (exclusive)");
        exit(2);
    }

    if let Some(t) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .expect("failed to configure thread pool");
    }

    let size_tag = if args.population_size {
        "popsize".to_string()
    } else {
        format!("n{}", args.size.unwrap())
    };
    let size_desc = if args.population_size {
        "population size".to_string()
    } else {
        format!("size {}", args.size.unwrap())
    };
    // Column tag for the CI bounds, e.g. 95 → "ci95", 99 → "ci99".
    let ci_tag = format!("ci{}", fmt_ci_level(args.ci_level));

    if args.carrier_mad {
        let groups = match load_carrier_groups(&args.input, &args.loci) {
            Ok(g) => g,
            Err(e) => {
                eprintln!("error: {e}");
                exit(1);
            }
        };
        if groups.is_empty() {
            report_empty_groups(&args.loci);
        }
        run_carrier_mad(&args, &groups, &size_tag, &size_desc, &ci_tag);
        return;
    }

    let groups = match load_groups(&args.input, &args.loci, args.allele_ci) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("error: {e}");
            exit(1);
        }
    };

    if groups.is_empty() {
        report_empty_groups(&args.loci);
    }

    if args.allele_ci {
        run_allele_ci(&args, &groups, &size_tag, &size_desc, &ci_tag);
    } else {
        run_mad(&args, &groups, &size_tag, &size_desc, &ci_tag);
    }
}

/// Print the "no groups found" error and exit(1); shared by every mode's
/// empty-groups check.
fn report_empty_groups(loci: &[String]) -> ! {
    if loci.is_empty() {
        eprintln!("No (level, code, locus) groups with any observations found.");
    } else {
        eprintln!(
            "No (level, code, locus) groups found for the requested loci: {:?}",
            loci
        );
    }
    exit(1);
}

/// MAD analysis: resample, compute MAD per cutoff, and summarise per
/// (level, code, locus, cutoff).
fn run_mad(args: &Args, groups: &[Group], size_tag: &str, size_desc: &str, ci_tag: &str) {
    let cutoffs = cutoffs();

    // Fan out over (rep, group) tasks. Each task is independent and deterministic.
    let mut rows: Vec<Row> = (0..args.reps)
        .into_par_iter()
        .flat_map_iter(|rep| {
            simulate_rep(rep, groups, &cutoffs, args.size, args.population_size, args.seed)
        })
        .collect();

    // Deterministic, stable ordering for output (rep, level, code, locus, cutoff).
    rows.par_sort_by(|a, b| {
        a.rep
            .cmp(&b.rep)
            .then_with(|| a.level.cmp(b.level))
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
            .then_with(|| a.cutoff.partial_cmp(&b.cutoff).unwrap())
    });

    let out = args
        .output
        .clone()
        .unwrap_or_else(|| format!("mad_simulations_{size_tag}.csv"));
    let summary_out = replace_csv_suffix(&out, "_summary.csv");

    if let Err(e) = write_sim_csv(&out, &rows) {
        eprintln!("error writing {out}: {e}");
        exit(1);
    }

    let summary = summarize_generic(
        &rows,
        &cutoffs,
        args.ci_level,
        |r| r.level,
        |r| r.code.as_str(),
        |r| r.locus.as_str(),
        |r| r.cutoff,
        |r| r.n_alleles,
    );
    if let Err(e) = write_summary_csv(&summary_out, &summary, ci_tag) {
        eprintln!("error writing {summary_out}: {e}");
        exit(1);
    }

    println!(
        "\nWrote {out} ({} rows; {} reps, {size_desc})",
        rows.len(),
        args.reps
    );
    println!("Wrote {summary_out}");
    println!("\nMean MAD per (level, code, locus, cutoff) [head]:");
    print_summary_head(&summary, 20, ci_tag);
}

/// Carrier-frequency MAD analysis: bootstrap whole samples (not individual
/// allele draws), then run the greedy carrier-coverage set cover to find the
/// minimum number of alleles needed to cover each cutoff fraction of people.
fn run_carrier_mad(
    args: &Args,
    groups: &[CarrierGroup],
    size_tag: &str,
    size_desc: &str,
    ci_tag: &str,
) {
    let cutoffs = cutoffs();

    // Fan out over (rep, group) tasks. Each task is independent and deterministic.
    let mut rows: Vec<CarrierRow> = (0..args.reps)
        .into_par_iter()
        .flat_map_iter(|rep| {
            simulate_rep_carrier(rep, groups, &cutoffs, args.size, args.population_size, args.seed)
        })
        .collect();

    // Deterministic, stable ordering for output (rep, level, code, locus, cutoff).
    rows.par_sort_by(|a, b| {
        a.rep
            .cmp(&b.rep)
            .then_with(|| a.level.cmp(b.level))
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
            .then_with(|| a.cutoff.partial_cmp(&b.cutoff).unwrap())
    });

    let out = args
        .output
        .clone()
        .unwrap_or_else(|| format!("carrier_mad_simulations_{size_tag}.csv"));
    let summary_out = replace_csv_suffix(&out, "_summary.csv");

    if let Err(e) = write_carrier_sim_csv(&out, &rows) {
        eprintln!("error writing {out}: {e}");
        exit(1);
    }

    let summary = summarize_generic(
        &rows,
        &cutoffs,
        args.ci_level,
        |r| r.level,
        |r| r.code.as_str(),
        |r| r.locus.as_str(),
        |r| r.cutoff,
        |r| r.n_alleles,
    );
    if let Err(e) = write_summary_csv(&summary_out, &summary, ci_tag) {
        eprintln!("error writing {summary_out}: {e}");
        exit(1);
    }

    println!(
        "\nWrote {out} ({} rows; {} reps, {size_desc})",
        rows.len(),
        args.reps
    );
    println!("Wrote {summary_out}");
    println!("\nMean carrier MAD per (level, code, locus, cutoff) [head]:");
    print_summary_head(&summary, 20, ci_tag);
}

/// Allele-frequency CI analysis: resample, record the bootstrap frequency of
/// every observed allele, and summarise per (level, code, locus, haplotype).
fn run_allele_ci(args: &Args, groups: &[Group], size_tag: &str, size_desc: &str, ci_tag: &str) {
    // Fan out over (rep, group) tasks. Each task is independent and deterministic.
    let mut rows: Vec<AlleleRow> = (0..args.reps)
        .into_par_iter()
        .flat_map_iter(|rep| {
            simulate_rep_allele_ci(rep, groups, args.size, args.population_size, args.seed)
        })
        .collect();

    // Deterministic, stable ordering for output (rep, level, code, locus, haplotype).
    rows.par_sort_by(|a, b| {
        a.rep
            .cmp(&b.rep)
            .then_with(|| a.level.cmp(b.level))
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
            .then_with(|| a.haplotype.cmp(&b.haplotype))
    });

    let out = args
        .output
        .clone()
        .unwrap_or_else(|| format!("allele_ci_simulations_{size_tag}.csv"));
    let summary_out = replace_csv_suffix(&out, "_summary.csv");

    if let Err(e) = write_allele_sim_csv(&out, &rows) {
        eprintln!("error writing {out}: {e}");
        exit(1);
    }

    let summary = summarize_allele_ci(&rows, args.ci_level);
    if let Err(e) = write_allele_summary_csv(&summary_out, &summary, ci_tag) {
        eprintln!("error writing {summary_out}: {e}");
        exit(1);
    }

    println!(
        "\nWrote {out} ({} rows; {} reps, {size_desc})",
        rows.len(),
        args.reps
    );
    println!("Wrote {summary_out}");
    println!("\nBootstrap allele-frequency CI per (level, code, locus, haplotype) [head]:");
    print_allele_summary_head(&summary, 20, ci_tag);
}

/// Resample size for a group: population-size mode uses the group's own observed
/// count; fixed mode uses `--size`, except the pooled `global` level which uses
/// `--size` x (number of pop2 populations for its locus).
fn resample_size(g: &Group, fixed_size: Option<usize>, population_size: bool) -> usize {
    if population_size {
        g.obs.len()
    } else if g.level == GLOBAL_LEVEL {
        fixed_size.unwrap() * g.n_pops
    } else {
        fixed_size.unwrap()
    }
}

/// Run every group for a single replicate, returning its output rows.
fn simulate_rep(
    rep: u64,
    groups: &[Group],
    cutoffs: &[f64],
    fixed_size: Option<usize>,
    population_size: bool,
    seed: u64,
) -> Vec<Row> {
    let mut out = Vec::with_capacity(groups.len() * cutoffs.len());
    for (g_idx, g) in groups.iter().enumerate() {
        if g.obs.is_empty() {
            continue;
        }
        let size = resample_size(g, fixed_size, population_size);

        // Deterministic per-(rep, group) seed → reproducible regardless of threading.
        let mut rng = ChaCha8Rng::seed_from_u64(mix_seed(seed, rep, g_idx as u64));

        // Sample with replacement, tally per-allele counts on the fly.
        let mut counts = vec![0u32; g.n_alleles];
        for _ in 0..size {
            let pick = g.obs[rng.gen_range(0..g.obs.len())] as usize;
            counts[pick] += 1;
        }

        // Allele frequencies of observed alleles, sorted descending.
        let mut pop_af: Vec<f64> = counts
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| c as f64 / size as f64)
            .collect();
        pop_af.sort_unstable_by(|a, b| b.partial_cmp(a).unwrap());
        let n_sampled_alleles = pop_af.len();

        // Cumulative AF, shared across all cutoffs for this sample.
        let mut cum: Vec<f64> = Vec::with_capacity(pop_af.len());
        let mut acc = 0.0;
        for &af in &pop_af {
            acc += af;
            cum.push(acc);
        }

        for &cutoff in cutoffs {
            let (n_alleles, cum_af) = mad_for_cutoff(&cum, cutoff);
            out.push(Row {
                rep,
                level: g.level,
                code: g.code.clone(),
                locus: g.locus.clone(),
                cutoff,
                n_alleles,
                cumulative_pop_af: round_half_even(cum_af, 6),
                n_sampled_alleles,
                sample_size: size,
            });
        }
    }
    out
}

/// Run every group for a single replicate in allele-CI mode, recording the
/// resampled frequency of *every* observed allele (carrying zeros for alleles
/// not drawn this rep so the CI reflects their full bootstrap distribution).
fn simulate_rep_allele_ci(
    rep: u64,
    groups: &[Group],
    fixed_size: Option<usize>,
    population_size: bool,
    seed: u64,
) -> Vec<AlleleRow> {
    let mut out = Vec::new();
    for (g_idx, g) in groups.iter().enumerate() {
        if g.obs.is_empty() {
            continue;
        }
        let size = resample_size(g, fixed_size, population_size);

        // Deterministic per-(rep, group) seed → reproducible regardless of threading.
        // Identical to the MAD path so both analyses see the same resamples.
        let mut rng = ChaCha8Rng::seed_from_u64(mix_seed(seed, rep, g_idx as u64));

        let mut counts = vec![0u32; g.n_alleles];
        for _ in 0..size {
            let pick = g.obs[rng.gen_range(0..g.obs.len())] as usize;
            counts[pick] += 1;
        }

        for id in 0..g.n_alleles {
            out.push(AlleleRow {
                rep,
                level: g.level,
                code: g.code.clone(),
                locus: g.locus.clone(),
                haplotype: g.labels[id].clone(),
                allele_af: round_half_even(counts[id] as f64 / size as f64, 6),
                sample_size: size,
            });
        }
    }
    out
}

/// Run every group for a single replicate in carrier-MAD mode: bootstrap whole
/// people with replacement, then greedily cover them by carried allele.
fn simulate_rep_carrier(
    rep: u64,
    groups: &[CarrierGroup],
    cutoffs: &[f64],
    fixed_size: Option<usize>,
    population_size: bool,
    seed: u64,
) -> Vec<CarrierRow> {
    let mut out = Vec::with_capacity(groups.len() * cutoffs.len());
    for (g_idx, g) in groups.iter().enumerate() {
        if g.persons.is_empty() {
            continue;
        }
        let size = if population_size {
            g.persons.len()
        } else {
            fixed_size.unwrap()
        };

        // Deterministic per-(rep, group) seed → reproducible regardless of threading.
        let mut rng = ChaCha8Rng::seed_from_u64(mix_seed(seed, rep, g_idx as u64));

        // Bootstrap draw of people; each draw carries that person's full
        // genotype at this locus into the resample.
        let draws: Vec<usize> = (0..size)
            .map(|_| rng.gen_range(0..g.persons.len()))
            .collect();

        let (steps, n_available_alleles) =
            greedy_carrier_cover(&draws, &g.persons, g.n_alleles, size);

        for &cutoff in cutoffs {
            let (n_alleles, frac) = carrier_mad_for_cutoff(&steps, cutoff);
            out.push(CarrierRow {
                rep,
                level: g.level,
                code: g.code.clone(),
                locus: g.locus.clone(),
                cutoff,
                n_alleles,
                cumulative_carrier_frac: round_half_even(frac, 6),
                n_available_alleles,
                sample_size: size,
            });
        }
    }
    out
}

/// Greedy set cover over a bootstrap draw of people: repeatedly pick the
/// allele carried by the most not-yet-covered people, remove every drawn
/// instance that carries it (their whole genotype, not just that allele), and
/// repeat until everyone is covered. Ports `carrier_mad.py`'s
/// `min_alleles_for_coverage`. Ties are broken by smallest allele id, for
/// determinism (the Python version's tie-break follows pandas groupby order,
/// which isn't a meaningful distinction to preserve).
///
/// Returns the cumulative-coverage step sequence `(alleles_chosen,
/// covered_frac)` after each pick — running to full coverage in one pass covers
/// every cutoff — plus the number of distinct alleles carried by the draw.
fn greedy_carrier_cover(
    draws: &[usize],
    persons: &[Vec<u32>],
    n_alleles: usize,
    size: usize,
) -> (Vec<(usize, f64)>, usize) {
    let mut active: Vec<usize> = (0..draws.len()).collect();
    let mut counts = vec![0u32; n_alleles];

    let mut available = std::collections::HashSet::new();
    for &idx in &active {
        for &a in &persons[draws[idx]] {
            available.insert(a);
        }
    }
    let n_available_alleles = available.len();

    let mut steps = Vec::new();
    let mut covered = 0usize;
    while !active.is_empty() {
        // Tally carrier counts over currently-active (not-yet-covered) draws.
        let mut touched: Vec<u32> = Vec::new();
        for &idx in &active {
            for &a in &persons[draws[idx]] {
                if counts[a as usize] == 0 {
                    touched.push(a);
                }
                counts[a as usize] += 1;
            }
        }

        let mut best_allele = touched[0];
        let mut best_count = counts[best_allele as usize];
        for &a in &touched {
            let c = counts[a as usize];
            if c > best_count || (c == best_count && a < best_allele) {
                best_count = c;
                best_allele = a;
            }
        }
        for &a in &touched {
            counts[a as usize] = 0;
        }

        let before = active.len();
        active.retain(|&idx| !persons[draws[idx]].contains(&best_allele));
        covered += before - active.len();
        steps.push((steps.len() + 1, covered as f64 / size as f64));
    }
    (steps, n_available_alleles)
}

/// Minimum number of alleles whose cumulative carrier coverage first exceeds
/// `cutoff`, from a `greedy_carrier_cover` step sequence.
fn carrier_mad_for_cutoff(steps: &[(usize, f64)], cutoff: f64) -> (usize, f64) {
    match steps.iter().position(|&(_, frac)| frac > cutoff) {
        Some(idx) => steps[idx],
        None => *steps.last().unwrap(),
    }
}

/// Minimum number of alleles whose cumulative AF first exceeds `cutoff`.
///
/// `cum` is the cumulative sum of the descending-sorted allele frequencies.
fn mad_for_cutoff(cum: &[f64], cutoff: f64) -> (usize, f64) {
    let total = cum.len();
    match cum.iter().position(|&c| c > cutoff) {
        Some(idx) => (idx + 1, cum[idx]),
        None => (total, *cum.last().unwrap()),
    }
}

/// Open a file as a byte reader, transparently decompressing gzip input.
///
/// Detection is by the gzip magic number (`0x1f 0x8b`), so `.csv.gz`, `.gz`,
/// or an extensionless gzip stream all work regardless of file name.
fn open_maybe_gzip(path: &str) -> Result<Box<dyn Read>, String> {
    let mut file = BufReader::new(
        File::open(path).map_err(|e| format!("could not open {path}: {e}"))?,
    );
    let mut magic = [0u8; 2];
    let n = file
        .read(&mut magic)
        .map_err(|e| format!("could not read {path}: {e}"))?;
    // Re-prepend the bytes we peeked so the chosen decoder sees a full stream.
    let chained = std::io::Cursor::new(magic[..n].to_vec()).chain(file);
    if n == 2 && magic == [0x1f, 0x8b] {
        Ok(Box::new(MultiGzDecoder::new(chained)))
    } else {
        Ok(Box::new(chained))
    }
}

/// Read all input CSVs and build the per-(level, code, locus) groups, ordered
/// by (level, code, locus) for deterministic output.
///
/// When `loci` is non-empty, only rows whose `locus` is in that set are kept.
/// When `allele_ci` is set, an extra `global` group per locus is built that
/// pools every population together, plus the per-locus pop2 population count
/// used to size its `--size`-mode resample.
fn load_groups(paths: &[String], loci: &[String], allele_ci: bool) -> Result<Vec<Group>, String> {
    // Optional locus allow-list; empty means "all loci".
    let locus_filter: Option<std::collections::HashSet<&str>> = if loci.is_empty() {
        None
    } else {
        Some(loci.iter().map(|s| s.as_str()).collect())
    };

    // Buckets keyed by (level_idx, code, locus). Value maps haplotype label →
    // dense allele id plus the running list of observations.
    struct Bucket {
        ids: HashMap<String, u32>,
        obs: Vec<u32>,
    }
    let mut buckets: HashMap<(usize, String, String), Bucket> = HashMap::new();

    // Sentinel level index for the pooled `global` group (`--allele-ci` only),
    // sitting just past the real per-population levels.
    let global_idx = LEVELS.len();
    // locus -> set of observed pop2 codes, for the global `--size` resample count.
    let mut pop2_by_locus: HashMap<String, std::collections::HashSet<String>> = HashMap::new();

    for path in paths {
        let inner = open_maybe_gzip(path)?;
        let mut reader = csv::Reader::from_reader(inner);
        let headers = reader
            .headers()
            .map_err(|e| format!("could not read header of {path}: {e}"))?
            .clone();

        let col = |name: &str| headers.iter().position(|h| h == name);
        let hap_i = col("haplotype");
        let locus_i = col("locus");
        let level_is: Vec<Option<usize>> = LEVELS.iter().map(|l| col(l)).collect();

        // Validate required columns (matches the Python missing-columns check).
        let mut missing: Vec<&str> = Vec::new();
        if hap_i.is_none() {
            missing.push("haplotype");
        }
        if locus_i.is_none() {
            missing.push("locus");
        }
        for (li, l) in LEVELS.iter().enumerate() {
            if level_is[li].is_none() {
                missing.push(l);
            }
        }
        if !missing.is_empty() {
            missing.sort_unstable();
            return Err(format!("Input is missing required columns: {missing:?}"));
        }
        let hap_i = hap_i.unwrap();
        let locus_i = locus_i.unwrap();

        let mut record = csv::StringRecord::new();
        while reader
            .read_record(&mut record)
            .map_err(|e| format!("error reading {path}: {e}"))?
        {
            let haplotype = record.get(hap_i).unwrap_or("");
            let locus = record.get(locus_i).unwrap_or("");
            // Skip loci outside the requested subset, if one was given.
            if let Some(filter) = &locus_filter {
                if !filter.contains(locus) {
                    continue;
                }
            }
            let locus = locus.to_string();
            for (level_idx, _level) in LEVELS.iter().enumerate() {
                let code = record.get(level_is[level_idx].unwrap()).unwrap_or("").to_string();
                let bucket = buckets
                    .entry((level_idx, code, locus.clone()))
                    .or_insert_with(|| Bucket {
                        ids: HashMap::new(),
                        obs: Vec::new(),
                    });
                let next_id = bucket.ids.len() as u32;
                let id = *bucket.ids.entry(haplotype.to_string()).or_insert(next_id);
                bucket.obs.push(id);
            }
            if allele_ci {
                // Track distinct pop2 codes per locus for the global size, and
                // pool every row into a single global bucket for the locus.
                let pop2_code = record.get(level_is[1].unwrap()).unwrap_or("").to_string();
                pop2_by_locus.entry(locus.clone()).or_default().insert(pop2_code);
                let bucket = buckets
                    .entry((global_idx, GLOBAL_LEVEL.to_string(), locus.clone()))
                    .or_insert_with(|| Bucket {
                        ids: HashMap::new(),
                        obs: Vec::new(),
                    });
                let next_id = bucket.ids.len() as u32;
                let id = *bucket.ids.entry(haplotype.to_string()).or_insert(next_id);
                bucket.obs.push(id);
            }
        }
    }

    // Distinct pop2 population count per locus, for the global `--size` resample.
    let n_pops_by_locus: HashMap<String, usize> = pop2_by_locus
        .into_iter()
        .map(|(locus, codes)| (locus, codes.len()))
        .collect();

    let mut groups: Vec<Group> = buckets
        .into_iter()
        .map(|((level_idx, code, locus), b)| {
            // Recover haplotype labels in dense-id order for allele-CI output.
            let mut labels = vec![String::new(); b.ids.len()];
            for (label, id) in &b.ids {
                labels[*id as usize] = label.clone();
            }
            // The sentinel index maps to the pooled `global` level; it also
            // carries the per-locus pop2 count for its `--size` resample.
            let (level, n_pops) = if level_idx < LEVELS.len() {
                (LEVELS[level_idx], 0)
            } else {
                (GLOBAL_LEVEL, *n_pops_by_locus.get(&locus).unwrap_or(&0))
            };
            Group {
                level,
                code,
                locus,
                n_alleles: b.ids.len(),
                obs: b.obs,
                labels,
                n_pops,
            }
        })
        .collect();

    groups.sort_by(|a, b| {
        a.level
            .cmp(b.level)
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
    });

    Ok(groups)
}

/// Read all input CSVs and build the per-(level, code, locus) carrier groups
/// for `--carrier-mad`, ordered by (level, code, locus) for deterministic
/// output. Unlike `load_groups`, rows are grouped by `sample` so each person's
/// full set of carried alleles at a locus stays together for the greedy
/// carrier-coverage cover; this requires an extra `sample` column beyond
/// `haplotype`/`locus`/`pop1`/`pop2`.
fn load_carrier_groups(paths: &[String], loci: &[String]) -> Result<Vec<CarrierGroup>, String> {
    // Optional locus allow-list; empty means "all loci".
    let locus_filter: Option<std::collections::HashSet<&str>> = if loci.is_empty() {
        None
    } else {
        Some(loci.iter().map(|s| s.as_str()).collect())
    };

    // Buckets keyed by (level_idx, code, locus). Value maps haplotype label →
    // dense allele id, plus each sample's set of carried allele ids.
    struct Bucket {
        ids: HashMap<String, u32>,
        persons: HashMap<String, std::collections::HashSet<u32>>,
    }
    let mut buckets: HashMap<(usize, String, String), Bucket> = HashMap::new();

    for path in paths {
        let inner = open_maybe_gzip(path)?;
        let mut reader = csv::Reader::from_reader(inner);
        let headers = reader
            .headers()
            .map_err(|e| format!("could not read header of {path}: {e}"))?
            .clone();

        let col = |name: &str| headers.iter().position(|h| h == name);
        let hap_i = col("haplotype");
        let locus_i = col("locus");
        let sample_i = col("sample");
        let level_is: Vec<Option<usize>> = LEVELS.iter().map(|l| col(l)).collect();

        // Validate required columns (matches the Python missing-columns check,
        // plus the extra `sample` column carrier mode needs).
        let mut missing: Vec<&str> = Vec::new();
        if hap_i.is_none() {
            missing.push("haplotype");
        }
        if locus_i.is_none() {
            missing.push("locus");
        }
        if sample_i.is_none() {
            missing.push("sample");
        }
        for (li, l) in LEVELS.iter().enumerate() {
            if level_is[li].is_none() {
                missing.push(l);
            }
        }
        if !missing.is_empty() {
            missing.sort_unstable();
            return Err(format!("Input is missing required columns: {missing:?}"));
        }
        let hap_i = hap_i.unwrap();
        let locus_i = locus_i.unwrap();
        let sample_i = sample_i.unwrap();

        let mut record = csv::StringRecord::new();
        while reader
            .read_record(&mut record)
            .map_err(|e| format!("error reading {path}: {e}"))?
        {
            let haplotype = record.get(hap_i).unwrap_or("");
            let locus = record.get(locus_i).unwrap_or("");
            // Skip loci outside the requested subset, if one was given.
            if let Some(filter) = &locus_filter {
                if !filter.contains(locus) {
                    continue;
                }
            }
            let locus = locus.to_string();
            let sample = record.get(sample_i).unwrap_or("").to_string();
            for (level_idx, _level) in LEVELS.iter().enumerate() {
                let code = record.get(level_is[level_idx].unwrap()).unwrap_or("").to_string();
                let bucket = buckets
                    .entry((level_idx, code, locus.clone()))
                    .or_insert_with(|| Bucket {
                        ids: HashMap::new(),
                        persons: HashMap::new(),
                    });
                let next_id = bucket.ids.len() as u32;
                let id = *bucket.ids.entry(haplotype.to_string()).or_insert(next_id);
                bucket.persons.entry(sample.clone()).or_default().insert(id);
            }
        }
    }

    let mut groups: Vec<CarrierGroup> = buckets
        .into_iter()
        .map(|((level_idx, code, locus), b)| {
            let persons: Vec<Vec<u32>> = b
                .persons
                .into_values()
                .map(|set| {
                    let mut v: Vec<u32> = set.into_iter().collect();
                    v.sort_unstable();
                    v
                })
                .collect();
            CarrierGroup {
                level: LEVELS[level_idx],
                code,
                locus,
                n_alleles: b.ids.len(),
                persons,
            }
        })
        .collect();

    groups.sort_by(|a, b| {
        a.level
            .cmp(b.level)
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
    });

    Ok(groups)
}

/// Aggregated summary stat per (level, code, locus, cutoff).
struct SummaryRow {
    level: &'static str,
    code: String,
    locus: String,
    cutoff: f64,
    mean: f64,
    median: f64,
    std: f64, // sample std (ddof=1); NaN when n < 2
    min: f64,
    max: f64,
    ci_low: f64,
    ci_high: f64,
}

/// Group n_alleles by (level, code, locus, cutoff) and compute summary stats,
/// mirroring the pandas groupby().agg() in the Python version. The reported CI
/// bounds are the two-sided percentiles for `ci_level` (e.g. 95 → 2.5/97.5).
///
/// Generic over the row type so both the allele-frequency `Row` and the
/// carrier-frequency `CarrierRow` (same per-cutoff `n_alleles` shape) share one
/// implementation.
fn summarize_generic<T>(
    rows: &[T],
    cutoffs: &[f64],
    ci_level: f64,
    level_of: impl Fn(&T) -> &'static str,
    code_of: impl Fn(&T) -> &str,
    locus_of: impl Fn(&T) -> &str,
    cutoff_of: impl Fn(&T) -> f64,
    n_alleles_of: impl Fn(&T) -> usize,
) -> Vec<SummaryRow> {
    let alpha = (100.0 - ci_level) / 2.0;
    let (low_pct, high_pct) = (alpha, 100.0 - alpha);
    // Key → collected n_alleles values.
    let mut groups: HashMap<(&'static str, String, String, u64), Vec<f64>> = HashMap::new();
    for r in rows {
        let key = (
            level_of(r),
            code_of(r).to_string(),
            locus_of(r).to_string(),
            cutoff_key(cutoff_of(r)),
        );
        groups.entry(key).or_default().push(n_alleles_of(r) as f64);
    }

    // Map back from quantized cutoff key to its float value.
    let cutoff_of: HashMap<u64, f64> = cutoffs.iter().map(|&c| (cutoff_key(c), c)).collect();

    let mut out: Vec<SummaryRow> = groups
        .into_iter()
        .map(|((level, code, locus, ck), mut vals)| {
            vals.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            let n = vals.len();
            let mean = vals.iter().sum::<f64>() / n as f64;
            let std = if n < 2 {
                f64::NAN
            } else {
                let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
                var.sqrt()
            };
            SummaryRow {
                level,
                code,
                locus,
                cutoff: cutoff_of[&ck],
                mean: round_half_even(mean, 3),
                median: round_half_even(percentile(&vals, 50.0), 3),
                std: round_half_even(std, 3),
                min: round_half_even(*vals.first().unwrap(), 3),
                max: round_half_even(*vals.last().unwrap(), 3),
                ci_low: round_half_even(percentile(&vals, low_pct), 3),
                ci_high: round_half_even(percentile(&vals, high_pct), 3),
            }
        })
        .collect();

    out.sort_by(|a, b| {
        a.level
            .cmp(b.level)
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
            .then_with(|| a.cutoff.partial_cmp(&b.cutoff).unwrap())
    });
    out
}

/// Aggregated allele-frequency summary per (level, code, locus, haplotype).
struct AlleleSummaryRow {
    level: &'static str,
    code: String,
    locus: String,
    haplotype: String,
    mean: f64,
    median: f64,
    std: f64, // sample std (ddof=1); NaN when n < 2
    min: f64,
    max: f64,
    ci_low: f64,
    ci_high: f64,
}

/// Group allele_af by (level, code, locus, haplotype) and compute summary stats,
/// mirroring the pandas groupby().agg() in the Python `--allele-ci` path. The
/// reported CI bounds are the two-sided percentiles for `ci_level`.
fn summarize_allele_ci(rows: &[AlleleRow], ci_level: f64) -> Vec<AlleleSummaryRow> {
    let alpha = (100.0 - ci_level) / 2.0;
    let (low_pct, high_pct) = (alpha, 100.0 - alpha);
    // Key → collected allele_af values.
    let mut groups: HashMap<(&'static str, String, String, String), Vec<f64>> = HashMap::new();
    for r in rows {
        let key = (r.level, r.code.clone(), r.locus.clone(), r.haplotype.clone());
        groups.entry(key).or_default().push(r.allele_af);
    }

    let mut out: Vec<AlleleSummaryRow> = groups
        .into_iter()
        .map(|((level, code, locus, haplotype), mut vals)| {
            vals.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
            let n = vals.len();
            let mean = vals.iter().sum::<f64>() / n as f64;
            let std = if n < 2 {
                f64::NAN
            } else {
                let var = vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n as f64 - 1.0);
                var.sqrt()
            };
            AlleleSummaryRow {
                level,
                code,
                locus,
                haplotype,
                mean: round_half_even(mean, 6),
                median: round_half_even(percentile(&vals, 50.0), 6),
                std: round_half_even(std, 6),
                min: round_half_even(*vals.first().unwrap(), 6),
                max: round_half_even(*vals.last().unwrap(), 6),
                ci_low: round_half_even(percentile(&vals, low_pct), 6),
                ci_high: round_half_even(percentile(&vals, high_pct), 6),
            }
        })
        .collect();

    out.sort_by(|a, b| {
        a.level
            .cmp(b.level)
            .then_with(|| a.code.cmp(&b.code))
            .then_with(|| a.locus.cmp(&b.locus))
            .then_with(|| a.haplotype.cmp(&b.haplotype))
    });
    out
}

/// Linear-interpolated percentile (NumPy's default "linear" method).
/// `sorted` must be ascending and non-empty.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let rank = (p / 100.0) * (n as f64 - 1.0);
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        return sorted[lo];
    }
    let frac = rank - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}

fn write_sim_csv(path: &str, rows: &[Row]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "rep,level,code,locus,cutoff,n_alleles,cumulative_pop_af,n_sampled_alleles,sample_size"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{}",
            r.rep,
            r.level,
            csv_field(&r.code),
            csv_field(&r.locus),
            fmt_num(r.cutoff),
            r.n_alleles,
            fmt_num(r.cumulative_pop_af),
            r.n_sampled_alleles,
            r.sample_size
        )?;
    }
    w.flush()
}

fn write_carrier_sim_csv(path: &str, rows: &[CarrierRow]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "rep,level,code,locus,cutoff,n_alleles,cumulative_carrier_frac,n_available_alleles,sample_size"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{}",
            r.rep,
            r.level,
            csv_field(&r.code),
            csv_field(&r.locus),
            fmt_num(r.cutoff),
            r.n_alleles,
            fmt_num(r.cumulative_carrier_frac),
            r.n_available_alleles,
            r.sample_size
        )?;
    }
    w.flush()
}

fn write_summary_csv(path: &str, rows: &[SummaryRow], ci_tag: &str) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "level,code,locus,cutoff,mean,median,std,min,max,{ci_tag}_low,{ci_tag}_high"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{}",
            r.level,
            csv_field(&r.code),
            csv_field(&r.locus),
            fmt_num(r.cutoff),
            fmt_num(r.mean),
            fmt_num(r.median),
            fmt_nan(r.std),
            fmt_num(r.min),
            fmt_num(r.max),
            fmt_num(r.ci_low),
            fmt_num(r.ci_high),
        )?;
    }
    w.flush()
}

fn print_summary_head(rows: &[SummaryRow], n: usize, ci_tag: &str) {
    println!(
        "{:>5} {:>8} {:>8} {:>6} {:>8} {:>8} {:>8} {:>6} {:>6} {:>9} {:>10}",
        "level", "code", "locus", "cutoff", "mean", "median", "std", "min", "max",
        format!("{ci_tag}_low"),
        format!("{ci_tag}_high")
    );
    for r in rows.iter().take(n) {
        println!(
            "{:>5} {:>8} {:>8} {:>6} {:>8} {:>8} {:>8} {:>6} {:>6} {:>9} {:>10}",
            r.level,
            r.code,
            r.locus,
            fmt_num(r.cutoff),
            fmt_num(r.mean),
            fmt_num(r.median),
            fmt_nan(r.std),
            fmt_num(r.min),
            fmt_num(r.max),
            fmt_num(r.ci_low),
            fmt_num(r.ci_high),
        );
    }
}

fn write_allele_sim_csv(path: &str, rows: &[AlleleRow]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "rep,level,code,locus,haplotype,allele_af,sample_size")?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{}",
            r.rep,
            r.level,
            csv_field(&r.code),
            csv_field(&r.locus),
            csv_field(&r.haplotype),
            fmt_num(r.allele_af),
            r.sample_size
        )?;
    }
    w.flush()
}

fn write_allele_summary_csv(
    path: &str,
    rows: &[AlleleSummaryRow],
    ci_tag: &str,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "level,code,locus,haplotype,mean,median,std,min,max,{ci_tag}_low,{ci_tag}_high"
    )?;
    for r in rows {
        writeln!(
            w,
            "{},{},{},{},{},{},{},{},{},{},{}",
            r.level,
            csv_field(&r.code),
            csv_field(&r.locus),
            csv_field(&r.haplotype),
            fmt_num(r.mean),
            fmt_num(r.median),
            fmt_nan(r.std),
            fmt_num(r.min),
            fmt_num(r.max),
            fmt_num(r.ci_low),
            fmt_num(r.ci_high),
        )?;
    }
    w.flush()
}

fn print_allele_summary_head(rows: &[AlleleSummaryRow], n: usize, ci_tag: &str) {
    println!(
        "{:>5} {:>8} {:>8} {:>10} {:>10} {:>10} {:>10} {:>8} {:>8} {:>11} {:>12}",
        "level", "code", "locus", "haplotype", "mean", "median", "std", "min", "max",
        format!("{ci_tag}_low"),
        format!("{ci_tag}_high")
    );
    for r in rows.iter().take(n) {
        println!(
            "{:>5} {:>8} {:>8} {:>10} {:>10} {:>10} {:>10} {:>8} {:>8} {:>11} {:>12}",
            r.level,
            r.code,
            r.locus,
            r.haplotype,
            fmt_num(r.mean),
            fmt_num(r.median),
            fmt_nan(r.std),
            fmt_num(r.min),
            fmt_num(r.max),
            fmt_num(r.ci_low),
            fmt_num(r.ci_high),
        );
    }
}

/// Stable hash combine (splitmix64) for deterministic per-task seeds.
fn mix_seed(seed: u64, rep: u64, group: u64) -> u64 {
    let mut z = seed
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(rep.wrapping_mul(0xD1B5_4A32_D192_ED03))
        .wrapping_add(group.wrapping_mul(0xCBF2_9CE4_8422_5C53));
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Quantize a cutoff to an integer key (hundredths) for hashing/grouping.
fn cutoff_key(c: f64) -> u64 {
    (c * 100.0).round() as u64
}

/// Round half to even (banker's rounding), matching Python's `round` and
/// pandas' `.round()`.
fn round_half_even(x: f64, ndigits: i32) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let factor = 10f64.powi(ndigits);
    let scaled = x * factor;
    let floor = scaled.floor();
    let diff = scaled - floor;
    let rounded = if (diff - 0.5).abs() < 1e-9 {
        // Exactly halfway → round to even.
        if (floor as i64) % 2 == 0 {
            floor
        } else {
            floor + 1.0
        }
    } else {
        scaled.round()
    };
    rounded / factor
}

/// Format a number with shortest round-trippable representation (Rust default),
/// rendering whole floats without a trailing `.0`-style mismatch is left as-is.
fn fmt_num(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        // Integral value: still print as the float pandas would (e.g. "2.0").
        // But cutoffs / counts read better as compact decimals; pandas writes
        // 2.0 for float columns, so keep one decimal for whole floats.
        format!("{:.1}", x)
    } else {
        format!("{}", x)
    }
}

/// Render the confidence level for a column tag: whole numbers drop the
/// decimal (95.0 → "95"), fractional levels keep it (97.5 → "97.5").
fn fmt_ci_level(level: f64) -> String {
    if level == level.trunc() {
        format!("{}", level as i64)
    } else {
        format!("{}", level)
    }
}

fn fmt_nan(x: f64) -> String {
    if x.is_nan() {
        String::new()
    } else {
        fmt_num(x)
    }
}

/// Minimal CSV field quoting (handles commas, quotes, newlines).
fn csv_field(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Replace a trailing `.csv` with `suffix` (e.g. `_summary.csv`), matching the
/// Python `out.replace(".csv", ...)` behavior for the common case.
fn replace_csv_suffix(out: &str, suffix: &str) -> String {
    if let Some(stripped) = out.strip_suffix(".csv") {
        format!("{stripped}{suffix}")
    } else {
        format!("{out}{suffix}")
    }
}
