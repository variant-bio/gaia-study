# Global Atlas of Immune Alleles (GAIA) Study

This repository accompanies the GAIA study. It contains the Nextflow pipelines
used to genotype immune loci and benchmark the calls, and the Rust tool used for
the majority allelic diversity (MAD) resampling analyses.

| Path | Contents |
|------|----------|
| [pipelines/vb_IMM.nf](pipelines/vb_IMM.nf) | Genotyping of HLA, IGHC (IGHG), FCGR and IG/TR V genes from WGS/WES CRAMs |
| [pipelines/vb_HAP.nf](pipelines/vb_HAP.nf) | hap-eval: base-level benchmarking of called alleles against truth haplotypes |
| [pipelines/modules/vb_pipeline.nf](pipelines/modules/vb_pipeline.nf) | Process definitions shared by both workflows |
| [pipelines/bin/](pipelines/bin/) | Python helpers the processes call (added to `PATH` by Nextflow) |
| [pipelines/configs/](pipelines/configs/) | A local Nextflow config with placeholder container images |
| [pipelines/examples/](pipelines/examples/) | Example sample sheets, BED and params files for both workflows |
| [analysis_scripts/](analysis_scripts/) | `mad_simulations`, a multithreaded Rust bootstrap of MAD and allele-frequency CIs |
| [environment.yml](environment.yml) | Conda environment for Nextflow, the helper scripts and the Rust toolchain |

---

## 1. Environment

### Conda (recommended)

```sh
conda env create -f environment.yml
conda activate gaia
```

This installs:

* **Nextflow** (≥ 23.10) and Java 17+
* **samtools** and **Python 3.10+** with `pysam`, `numpy`, `pandas`, `scipy`,
  which is everything `vb_HAP.nf` and [pipelines/bin/](pipelines/bin/) need
* **Rust** (≥ 1.85) for building `mad_simulations`

### Without conda

| Need | Install |
|------|---------|
| Nextflow | `curl -s https://get.nextflow.io \| bash` (needs Java 17+) |
| Python deps | `pip install pysam numpy pandas scipy` |
| samtools | your package manager, e.g. `brew install samtools` / `apt install samtools` |
| Rust | [rustup](https://rustup.rs): `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \| sh` |

### Containers (vb_IMM.nf only)

The genotypers that `vb_IMM.nf` runs are not in the conda environment. Each
process expects a container image with its tool at a fixed path, because the
process scripts call those paths directly. To run `vb_IMM.nf` you'll need to
build images that match this layout (or edit the scripts in
[vb_pipeline.nf](pipelines/modules/vb_pipeline.nf)):

| Process | Tool | Paths and commands the script expects |
|---------|------|----------------------------------|
| `hlaTypeT1K`, `ighgTypeT1K` | [T1K](https://github.com/mourisl/T1K) 1.0.7 | `/root/T1K-1.0.7/run-t1k`; `aws` CLI on `PATH` |
| `hlaType` | [HLA*LA](https://github.com/DiltheyLab/HLA-LA) | `/root/HLA-LA/src/HLA-LA.pl`, `/root/picard.jar`, `/db/Homo_sapiens_assembly38.fasta` |
| `hlaTypeKourami` | [Kourami](https://github.com/Kingsford-Group/kourami) | `/root/Kourami.jar`, `$KOURAMI_DIR/scripts/alignAndExtract_hs38DH.sh`, `/db/Homo_sapiens_assembly38.fasta` |
| `immunotyper2` | [ImmunoTyper-SR](https://github.com/WGLab/ImmunoTyper-SR) | conda env `immunotyper-SR` at `/opt/conda`; `/db/Homo_sapiens_assembly38.fasta`; large scratch at `/tmp` |
| `locityperMain` | [Locityper](https://github.com/tprodanov/locityper) + jellyfish | `locityper` on `PATH` |
| `subsetBamFiles`, `compareCigarStrings`, `findMatches`, `summarizeResults` | hap-eval scripts | `python3` + `pysam`/`pandas`/`scipy`, `samtools` (the `gaia` conda env works as is) |

The CRAM processes also set `REF_PATH`/`REF_CACHE=/db/ref_cache/%2s/%2s/%s`,
so the images need a populated htslib reference cache for GRCh38 (built with
`seq_cache_populate.pl` from samtools). Without it, CRAM decoding falls back to
downloading reference sequences from the EBI.

---

## 2. Building the Rust binary (`mad_simulations`)

```sh
cd analysis_scripts
cargo build --release --locked
./target/release/mad_simulations --help
```

* `--locked` builds against the committed `Cargo.lock`, so you get exactly the
  dependency versions used for the paper.
* The release profile already enables LTO and `codegen-units = 1`. To tune for
  the build machine's CPU as well:
  `RUSTFLAGS="-C target-cpu=native" cargo build --release --locked`.
* To put the binary on your `PATH`:
  `cargo install --path analysis_scripts --locked` (installs to `~/.cargo/bin`).
* For a fully static Linux binary that runs on any cluster node:
  ```sh
  rustup target add x86_64-unknown-linux-musl
  cargo build --release --locked --target x86_64-unknown-linux-musl
  # -> target/x86_64-unknown-linux-musl/release/mad_simulations
  ```

Quick check:

```sh
printf 'haplotype,pop1,pop2,locus\nA*01,EUR,CEU,HLA-A\nA*02,EUR,CEU,HLA-A\nA*01,EUR,GBR,HLA-A\nA*03,AFR,YRI,HLA-A\n' > toy.csv
./target/release/mad_simulations toy.csv -n 10 -r 20 -s 1
# -> mad_simulations_n10.csv, mad_simulations_n10_summary.csv
```

[analysis_scripts/README.md](analysis_scripts/README.md) covers the input
format, the three analysis modes (MAD, `--allele-ci`, `--carrier-mad`), every
option, and the output columns.

---

## 3. Pipelines

All commands below run from [pipelines/](pipelines/). Sample-sheet paths must
be **absolute** paths or object-store URIs (`s3://`, `gs://`, …).

> **Set boolean parameters in the params file, not on the command line.**
> Recent Nextflow versions pass `--runHLA false` to the script as the string
> `"false"`, which is truthy, so the step runs anyway.

### Config

[configs/nextflow.local.config](pipelines/configs/nextflow.local.config) runs
every process with the local executor and sets per-process CPU and memory. It
has `docker` and `singularity` profiles. Its container images are
placeholders (`REPLACE_ME/...`): replace them with your own builds (see
[Containers](#containers-vb_immnf-only)). To run on a cluster, add the
executor settings your scheduler needs (e.g. `process.executor = 'slurm'`).

### 3a. vb_IMM.nf: immune-locus genotyping

Genotypes each sample from an aligned GRCh38 CRAM/BAM with the tools
enabled in the params file:

| Flag | Tool(s) | Output folder per sample |
|------|---------|--------------------------|
| `runHLA` + `hlaCallers` | T1K (`hla-wgs` preset, or `hla` if `exome: true`), HLA*LA, Kourami | `hlaT1K/`, `hlaOut/`, (Kourami: work dir only) |
| `runIGHG` | T1K with the IGHC database | `ighgT1K/` |
| `runFCGR` | Locityper against GRCh38 and CHM13 (needs FASTQs) | `locityper/GRCh38/`, `locityper/hs1/` |
| `runIGTRV` | ImmunoTyper-SR, one task per locus in `igtr_loci` | `immunotyper/` |

Results are published to
`<outputBucket>/<studyID>/<libraryTypeCode>/<sampleID>/`.

**Inputs** (examples in [pipelines/examples/imm/](pipelines/examples/imm/)):

* [`input.txt`](pipelines/examples/imm/input.txt): tab-separated, one row per sample
  ```
  sampleID	cramfile	craifile
  HG00096	/data/gaia/cram/HG00096.final.cram	/data/gaia/cram/HG00096.final.cram.crai
  ```
  BAMs also work (`.bam` + `.bai` in the same columns). CRAMs must be aligned
  to GRCh38. The 1000 Genomes high-coverage CRAMs are public at
  `s3://1000genomes/1000G_2504_high_coverage/` and on the EBI FTP.
* [`input_fastq.txt`](pipelines/examples/imm/input_fastq.txt): tab-separated, one row per
  read pair; rows sharing a `sampleID` are grouped (e.g. multiple lanes).
  `.gz` and `.bz2` both work. It's only used when `runFCGR: true`, but the file must
  exist (header-only is fine).
* [`params.yaml`](pipelines/examples/imm/params.yaml): run flags, output
  location and every resource path, with comments.

**Resources** (set in `params.yaml`):

* T1K allele databases: `hla_t1k_db.tar.gz` and `ighc_t1k_db.tar.gz` on
  Zenodo (Materials S1 and S2). The T1K processes copy them with `aws s3 cp` from
  `<resourceBucket>/pipeline/hla/t1k/hlaidx_dna_{seq,coord}.fa` and
  `<resourceBucket>/pipeline/t1k-ighg/HPRC-R2_IGHG_dna_{seq,coord_v2}.fa`, so
  upload them to that layout in your own bucket.
* HLA*LA graph `PRG_MHC_GRCh38_withIMGT` (`hlaDB`), Kourami panel (`kouramiDB`).
* Locityper: GRCh38 no-alt FASTA and CHM13 (`hs1`) FASTA with `.fai`, jellyfish
  k-mer counts and Locityper databases for each. The module also has
  `jellyfishRef` and `locityperAddDb` processes for building these once. Even
  when `runFCGR: false`, these params must point at existing files.

**Run:**

```sh
nextflow run vb_IMM.nf \
    -c configs/nextflow.local.config -profile docker \
    -params-file examples/imm/params.yaml
```

`studyID` and `libraryTypeCode` are required; the workflow exits with status 1
and no message if either is missing.

### 3b. vb_HAP.nf: hap-eval benchmarking

Scores called alleles against truth haplotypes base by base over each BED
region, as used for Tables S2, S20, S24, S25, S27 and S29:

1. `subsetBamFiles` (`trim_cigar.py`): trims every alignment in each BAM to
   each region
2. `compareCigarStrings` (`cigar_comp.py`): compares every called allele with
   every truth haplotype and counts TP/TN/FP/FN per aligned reference base
3. `findMatches` (`find_mates.py`): pairs called alleles with truth
   haplotypes using minimum-error Hungarian assignment, keeping only pairs with
   ≥ 500 observed bases
4. `summarizeResults` (`summarize_csv.py`): writes the summary tables
   `sum`, `per_hap`, `cov_reg`, `sample_cov`, `cov`, `sample`, `reg` (`.csv`) to
   `<outputBucket>/<studyID>/hap-eval/<libraryTypeCode>/`

**Inputs** (examples in [pipelines/examples/hap/](pipelines/examples/hap/)):

* [`comp_manifest.csv`](pipelines/examples/hap/comp_manifest.csv): the calls to evaluate
  ```
  sample,id,bam
  HG00731,HG00731_10X,/data/gaia/calls/HG00731_10X.calls.bam
  ```
* [`truth_manifest.csv`](pipelines/examples/hap/truth_manifest.csv): the truth haplotypes
  ```
  sample,hap,bam
  HG00731,hap1,/data/gaia/truth/HG00731.hap1.bam
  HG00731,hap2,/data/gaia/truth/HG00731.hap2.bam
  ```
* [`regions.bed`](pipelines/examples/hap/regions.bed): headerless BED, 4th column is
  the region name (example coordinates are illustrative GRCh38 gene spans)
* [`params.yaml`](pipelines/examples/hap/params.yaml)

Requirements that the scripts don't check (they produce wrong numbers rather than errors):

* **BAMs must use `=`/`X` CIGAR operations** (e.g. `minimap2 -a --eqx`).
  Mismatches are detected from `X` operations, so with plain `M` CIGARs
  every SNV is scored as a match.
* Each BAM needs a `.bai` next to it (`<bam>.bai`).
* Calls and truth are joined on `sample`, so the `sample` column must match
  between the two manifests.
* `id` must be `<SAMPLE>_<COVERAGE>` (e.g. `HG00731_10X`), where neither part
  contains `_` or `-`. `find_mates.py` parses `SAMPLE` and `SEQ_COV` from it.
* Region names can't contain `-`, and `hap` values can't contain `-` or `.`.
* Truth BAMs are typically haplotype-resolved assembly contigs aligned to the
  reference; call BAMs are the called allele sequences aligned the same way.

**Run:**

```sh
conda activate gaia     # no container needed
nextflow run vb_HAP.nf \
    -c configs/nextflow.local.config \
    -params-file examples/hap/params.yaml
```

### Helper scripts

The scripts in [pipelines/bin/](pipelines/bin/) also run standalone (`--help`
on each):

| Script | Purpose |
|--------|---------|
| `trim_cigar.py -b in.bam -d regions.bed -o out.bam` | Trim alignments to BED regions |
| `cigar_comp.py -t truth.bam -c calls.bam -o out.csv [-g genes.gff3 -n GENE]` | Per-base comparison; `-g/-n` restricts to exon/CDS of one gene (coding-region results, Table S25) |
| `find_mates.py -i *.csv -o out.csv` | Optimal call↔truth pairing |
| `summarize_csv.py -i *.csv -o outdir/` | Summary tables |
| `compile_gates.py {sample,merge,study}` | Collects per-sample T1K / HLA*LA / ImmunoTyper outputs into per-study tables |

---

## 4. Assumptions in the code

The pipeline code makes these assumptions, which you may need to work around:

| Where | What | Workaround |
|-------|------|------------|
| [vb_pipeline.nf](pipelines/modules/vb_pipeline.nf) `hlaTypeT1K`, `ighgTypeT1K`, `hlaTypeT1K_rna` | Fetch T1K databases with `aws s3 cp ${params.resourceBucket}/pipeline/...` inside the task, which needs AWS credentials and the `aws` CLI in the image | host the files at the same layout in your own bucket |
| [vb_pipeline.nf](pipelines/modules/vb_pipeline.nf) (all genotypers) | Absolute in-image paths (`/root/...`, `/db/...`, `/opt/conda`) | match the layout in [Containers](#containers-vb_immnf-only) |
| [vb_pipeline.nf](pipelines/modules/vb_pipeline.nf) `hlaType` | Derives the HLA*LA sample name as the 3rd `-`-separated field of `sampleID` (`A-B-HG00731-1` → `HG00731`); IDs without hyphens are used as is, and IDs with exactly one hyphen give an empty name | use IDs without hyphens |
| [compile_gates.py](pipelines/bin/compile_gates.py) | Expects `STDY-<code>` study names, sample IDs containing `-<code>-`, and `LBTP-*` library-type folders | name outputs to match |

`LBTP-WGS-5X` / `-10X` / `-30X` in Tables S15–S19 are library-type codes for
the downsampled coverages (5×, 10×, 30× WGS).
