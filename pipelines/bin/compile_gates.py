#!/usr/bin/env python3
"""Compile typing results: per sample, then per study, then across cohorts.

vb_IMM.nf publishes results per sample under

    <root>/<library type>/<sample>/<tool folder>/...

and they are compiled in three steps, each a mode of this script:

    sample   one sample's tool folder(s) to tables
             compile_gates.py sample --study STDY-ABC --sample SAMP-ABC-0001-1 \\
                 --lbtp LBTP-WGS-30X --tool immunotyper [--dir .] [--loci F] [--genes F]
             -> <sample>.<lbtp>.<tool>.<table>.csv.gz, <sample>.<lbtp>.<tool>.samples.csv

    merge    concatenate the tables of many samples (-> a study) or of many
             studies (-> all cohorts), table by table
             compile_gates.py merge --name STDY-ABC --tool immunotyper FILES... \\
                 [--manifest M --study STDY-ABC]
             -> <name>.<tool>.<table>.csv.gz, <name>.<tool>.samples.csv
             With --manifest, every row of the study that has no sample results
             (no tool folder, or its task failed) gets a "no results" status row.

    study    the whole of one study in one go, from a local results folder
             compile_gates.py study --root results/GAIA \\
                 --study ABC --tool immunotyper [--manifest M]

Every table is written on every run, header-only when empty, so the set of
output files for a tool never varies. samples.csv has one row per sample (per
locus, for tools with loci) saying what was found and what was missing.

Naming the code expects: study codes carry a STDY- prefix (added if missing),
sample IDs embed the study code as -<code>- (SAMP-ABC-0001-1 -> ABC), and
library-type folders start with LBTP-.

A manifest is a CSV with columns sample, lbtp and optionally study. A row's
study is its study column where filled in, else the code in its sample ID.

Tools (and the folders each reads):
    immunotyper   immunotyper/            (immunotyper2)
    t1k           hlaT1K/, ighgT1K/       (hlaTypeT1K, ighgTypeT1K)
    hla_la        hlaOut/                 (hlaType)
Adding one - locityper, from locityper/<ref>/<sample>/, say - means writing a
collector and adding it to COLLECTORS.
"""
import argparse
import csv
import gzip
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Optional, TextIO, cast

ID_COLUMNS = ["study", "sample", "lbtp"]


# --------------------------------------------------------------------------- #
# readers
# --------------------------------------------------------------------------- #
def data_lines(path, comment: Optional[str] = "#"):
    """The non-blank, non-comment lines of a text file (comment=None keeps all)."""
    with open(path) as fh:
        return [ln.rstrip("\n") for ln in fh
                if ln.strip() and not (comment and ln.startswith(comment))]


# On a network file system a read can fail transiently with EIO and succeed on
# a retry. Either way one file must not sink the rest, so a file that still
# fails is recorded against its sample and skipped.
READ_ATTEMPTS = 4
_tally = {"read": 0, "failed": 0}
_tally_lock = threading.Lock()


def read_lines(path, errors, comment: Optional[str] = "#"):
    """data_lines(), retried on I/O errors; None (and noted in `errors`) if the
    file cannot be read at all."""
    for attempt in range(READ_ATTEMPTS):
        try:
            lines = data_lines(path, comment)
            with _tally_lock:
                _tally["read"] += 1
            return lines
        except OSError as exc:
            if attempt + 1 < READ_ATTEMPTS:
                time.sleep(0.5 * 2 ** attempt)
                continue
            with _tally_lock:
                _tally["failed"] += 1
            errors.append(f"{os.path.basename(path)}: {exc.strerror or exc}")
            print(f"unreadable after {READ_ATTEMPTS} attempts: {path}: {exc}",
                  file=sys.stderr)
            return None


def names_in(path):
    """The names in a one-per-line file (or the first column of a CSV); a header
    line such as "gene" or "locus" is ignored."""
    if not path:
        return set()
    names = {ln.split(",")[0].strip() for ln in data_lines(path)}
    return names - {"gene", "genes", "locus", "loci", ""}


def fields(line, sep, n):
    """Split a line into exactly n fields (padding or folding the tail)."""
    parts = line.split(sep) if sep else line.split()
    if len(parts) > n:
        parts = parts[:n - 1] + [(sep or " ").join(parts[n - 1:])]
    return parts + [""] * (n - len(parts))


# --------------------------------------------------------------------------- #
# collectors: (folder holding the tool folders, sample ID, options)
#             -> ({table: [row, ...]}, [status row, ...])
# --------------------------------------------------------------------------- #
# immunotyper-SR, run once per locus into the same folder, writes
#   <cram basename>-<locus>_allele_calls.txt              one allele per line
#   <cram basename>-<locus>_functional_allele_calls.txt   one allele per line
#   <cram basename>-novel_variant_vcfs/<gene>_variants.vcf
# The basename is not fixed (downsampled CRAMs carry an extra suffix), so
# files are matched by suffix rather than rebuilt from the sample ID. The locus
# is written as given to --gene_type, which has been both ighv and TRBV, so it
# is matched in any case and reported in lower case.
IT_CALLS_RE = re.compile(
    r"^.+-(?P<locus>[a-z]+)_(?P<kind>functional_allele_calls|allele_calls)\.txt$",
    re.IGNORECASE)
VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
               "FORMAT", "SAMPLE"]


def gene_locus(gene):
    """IGHV1-2 -> ighv, TRAV12-1 -> trav."""
    return gene[:4].lower()


def collect_immunotyper(sdir, sample, opts):
    tables = {"allele_calls": [], "functional_allele_calls": [], "novel_variants": []}
    folder = os.path.join(sdir, "immunotyper")
    if not os.path.isdir(folder):
        return tables, [{"locus": "", "status": "no immunotyper folder"}]

    counts, errors = {}, {}

    def count(locus):
        return counts.setdefault(locus, {"allele_calls": 0, "functional_allele_calls": 0,
                                         "novel_genes": 0, "novel_variants": 0})

    unmatched = []
    for name in sorted(os.listdir(folder)):
        path = os.path.join(folder, name)
        m = IT_CALLS_RE.match(name)
        if not m and "allele_calls" in name.lower():
            unmatched.append(name)
        if m and os.path.isfile(path):
            locus, kind = m.group("locus").lower(), m.group("kind").lower()
            if opts.loci and locus not in opts.loci:
                continue
            alleles = read_lines(path, errors.setdefault(locus, []), comment=None)
            if alleles is None:
                continue
            tables[kind] += [{"locus": locus, "allele": a.strip()} for a in alleles]
            count(locus)[kind] += len(alleles)
        elif name.endswith("-novel_variant_vcfs") and os.path.isdir(path):
            for vcf in sorted(os.listdir(path)):
                if not vcf.endswith("_variants.vcf"):
                    continue
                gene = vcf[:-len("_variants.vcf")]
                locus = gene_locus(gene)
                if (opts.genes and gene not in opts.genes) or \
                        (opts.loci and locus not in opts.loci):
                    continue
                lines = read_lines(os.path.join(path, vcf), errors.setdefault(locus, []))
                if lines is None:
                    continue
                rows = [dict(zip(VCF_COLUMNS, fields(ln, "\t", len(VCF_COLUMNS))),
                             locus=locus, gene=gene)
                        for ln in lines]
                tables["novel_variants"] += rows
                count(locus)["novel_genes"] += 1
                count(locus)["novel_variants"] += len(rows)

    if unmatched:
        # looks like an allele-call file but is not named as expected: say so
        # rather than let it vanish from the tables
        print(f"{sample}: ignored {len(unmatched)} file(s) not named "
              f"<prefix>-<locus>_[functional_]allele_calls.txt: "
              f"{', '.join(unmatched[:5])}{' ...' if len(unmatched) > 5 else ''}",
              file=sys.stderr)

    # a locus task that failed leaves no allele calls: report the loci that were
    # expected and are absent, not only the ones that are there
    status = []
    for locus in sorted(set(counts) | set(opts.loci or []) | {k for k, v in errors.items() if v}):
        c = counts.get(locus)
        bad = errors.get(locus) or []
        if bad:
            state = "read error"
        elif not c or not c["allele_calls"]:
            state = "no allele calls"
        else:
            state = "ok"
        status.append(dict(c or count(locus), locus=locus, status=state,
                           unreadable="; ".join(bad)))
    return tables, status


T1K_GENOTYPE = ["locus", "cn", "allele_group_1", "abundance_1", "quality_1",
                "allele_group_2", "abundance_2", "quality_2", "secondary_alleles"]
T1K_ALLELE = ["allele", "quality"]
T1K_VCF = ["allele", "pos1", "dot1", "ref", "alt", "dot2", "qual", "alt_support",
           "read_depth", "uniq_support", "pos", "last"]


def collect_t1k(sdir, sample, opts):
    """HLA and IGHG, T1K: genotypes, alleles and the allele 'vcf'."""
    tables, status = {}, []
    for target in ("hla", "ighg"):
        folder = os.path.join(sdir, f"{target}T1K")
        found = {}
        for table, suffix, sep, columns in (
                ("genotypes", "genotype.tsv", "\t", T1K_GENOTYPE),
                ("alleles", "allele.tsv", None, T1K_ALLELE),
                ("allele_vcf", "allele.vcf", None, T1K_VCF)):
            path = os.path.join(folder, f"T1K_{sample}_{suffix}")
            rows = ([dict(zip(columns, fields(ln, sep, len(columns))))
                     for ln in data_lines(path)] if os.path.isfile(path) else [])
            tables[f"{target}_{table}"] = rows
            found[table] = len(rows) if os.path.isfile(path) else None
        state = ("no T1K folder" if not os.path.isdir(folder) else
                 "no genotype file" if found["genotypes"] is None else "ok")
        status.append(dict(locus=target, status=state,
                           **{k: ("" if v is None else v) for k, v in found.items()}))
    return tables, status


def collect_hla_la(sdir, sample, opts):
    path = os.path.join(sdir, "hlaOut", f"{sample}.R1_bestguess_G.txt")
    if not os.path.isfile(path):
        return {"bestguess_G": []}, [{"status": "no bestguess_G file", "calls": ""}]
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return {"bestguess_G": rows}, [{"status": "ok", "calls": len(rows)}]


# (collector, tables in output order with their columns - None = whatever the
# files have, in the order first seen - and the samples.csv columns)
COLLECTORS = {
    "immunotyper": (collect_immunotyper, {
        "allele_calls": ["locus", "allele"],
        "functional_allele_calls": ["locus", "allele"],
        "novel_variants": ["locus", "gene"] + VCF_COLUMNS,
    }, ["locus", "status", "allele_calls", "functional_allele_calls",
        "novel_genes", "novel_variants", "unreadable"]),
    "t1k": (collect_t1k, {
        f"{target}_{table}": columns
        for target in ("hla", "ighg")
        for table, columns in (("genotypes", T1K_GENOTYPE), ("alleles", T1K_ALLELE),
                               ("allele_vcf", T1K_VCF))
    }, ["locus", "status", "genotypes", "alleles", "allele_vcf", "unreadable"]),
    "hla_la": (collect_hla_la, {"bestguess_G": None}, ["status", "calls", "unreadable"]),
}


def compile_sample(folder, study, sample, lbtp, tool, opts):
    """One sample's tables and status rows, with the ID columns filled in."""
    collect = COLLECTORS[tool][0]
    ids = {"study": study, "sample": sample, "lbtp": lbtp}
    if not os.path.isdir(folder):          # a manifest row with no results
        return {}, [dict(ids, status="no sample folder")]
    try:
        tables, status = collect(folder, sample, opts)
    except OSError as exc:
        # listing a folder, or the T1K / HLA-LA readers: lose the sample only
        print(f"{sample}: {exc}", file=sys.stderr)
        with _tally_lock:
            _tally["failed"] += 1
        return {}, [dict(ids, status="read error", unreadable=str(exc))]
    return ({t: [dict(ids, **r) for r in rows] for t, rows in tables.items()},
            [dict(ids, **r) for r in status])


# --------------------------------------------------------------------------- #
# studies and manifests
# --------------------------------------------------------------------------- #
def as_study(code):
    code = code.strip()
    return code if code.startswith("STDY-") else f"STDY-{code}"


def study_code(study):
    """STDY-ABC -> ABC (the part sample IDs carry)."""
    return study[len("STDY-"):] if study.startswith("STDY-") else study


def row_study(row, column):
    """A manifest row's study code: its study column, else its sample ID."""
    given = (row.get(column) or "").strip()
    if given:
        return study_code(given)
    parts = row["sample"].split("-")
    return parts[1] if len(parts) > 2 else ""


def read_manifest(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = {"sample", "lbtp"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"{path}: needs columns sample and lbtp "
                     f"(has: {', '.join(reader.fieldnames or [])})")
        return list(reader)


def manifest_rows(path, code, column, library_types=None):
    """[(lbtp, sample)] for the manifest rows of one study, and every study named."""
    rows = read_manifest(path)
    out, seen = [], set()
    for row in rows:
        sample, lbtp = row["sample"].strip(), row["lbtp"].strip()
        if row_study(row, column) != code or (sample, lbtp) in seen:
            continue
        if library_types and lbtp not in library_types:
            continue
        seen.add((sample, lbtp))
        out.append((lbtp, sample))
    studies = sorted({row_study(r, column) for r in rows} - {""})
    return studies, out


def find_samples(root, code, library_types=None):
    """[(lbtp, sample)] for every sample folder of the study under root.

    Only LBTP-* folders are searched (that skips compiled/, hap-eval/, tmp/),
    and within them only names containing -<code>-.
    """
    tag = f"-{code}-"
    out = []
    with os.scandir(root) as top:
        lbtps = sorted(e.name for e in top if e.name.startswith("LBTP-") and e.is_dir())
    for lbtp in lbtps:
        if library_types and lbtp not in library_types:
            continue
        with os.scandir(os.path.join(root, lbtp)) as entries:
            for e in sorted(entries, key=lambda x: x.name):
                if tag in e.name and e.is_dir():
                    out.append((lbtp, e.name))
    return lbtps, out


# --------------------------------------------------------------------------- #
# writing and merging
# --------------------------------------------------------------------------- #
def open_text(path, mode) -> TextIO:
    """Text-mode open of a .csv or .csv.gz file."""
    if path.endswith(".gz"):
        return cast(TextIO, gzip.open(path, mode + "t", newline=""))
    return open(path, mode, newline="")


def write_csv(path, columns, rows):
    with open_text(path, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def union_columns(rows):
    seen = {}
    for row in rows:
        seen.update(dict.fromkeys(row))
    return list(seen)


def table_columns_for(tool, table, headers=()):
    """ID columns + the table's own, or for a free-form table whatever the
    inputs have."""
    fixed = COLLECTORS[tool][1].get(table)
    if fixed is not None:
        return ID_COLUMNS + fixed
    seen = dict.fromkeys(ID_COLUMNS)
    for header in headers:
        seen.update(dict.fromkeys(header))
    return list(seen)


def write_outputs(outdir, stem, tool, tables, status_rows):
    """<stem>.<tool>.<table>.csv.gz for every table the tool has, and
    <stem>.<tool>.samples.csv."""
    os.makedirs(outdir, exist_ok=True)
    for table in COLLECTORS[tool][1]:
        rows = tables.get(table, [])
        columns = table_columns_for(tool, table, [union_columns(rows)])
        path = os.path.join(outdir, f"{stem}.{tool}.{table}.csv.gz")
        write_csv(path, columns, rows)
        print(f"  {os.path.basename(path)}: {len(rows)} row(s)")
    path = os.path.join(outdir, f"{stem}.{tool}.samples.csv")
    write_csv(path, ID_COLUMNS + COLLECTORS[tool][2], status_rows)
    bad = sum(1 for r in status_rows if r.get("status") != "ok")
    print(f"  {os.path.basename(path)}: {len(status_rows)} row(s), {bad} not ok")
    return path


def table_of(path, tool):
    """<anything>.<tool>.<table>.csv[.gz] -> table ("samples" for samples.csv)."""
    name = os.path.basename(path)
    marker = f".{tool}."
    if marker not in name:
        return None
    rest = name[name.rindex(marker) + len(marker):]
    for ext in (".csv.gz", ".csv"):
        if rest.endswith(ext):
            return rest[:-len(ext)]
    return None


def check_read_errors():
    if _tally["failed"]:
        print(f"\n{_tally['failed']} file(s) could not be read ({_tally['read']} were) - "
              f"see the unreadable column of samples.csv.", file=sys.stderr)
        if not _tally["read"]:
            # nothing at all was readable: the tables are empty for a reason that
            # has nothing to do with the samples - fail rather than publish
            sys.exit("no file could be read - this is the storage or the "
                     "permissions, not the data")


# --------------------------------------------------------------------------- #
# modes
# --------------------------------------------------------------------------- #
def run_sample(opts):
    study = as_study(opts.study)
    tables, status = compile_sample(opts.dir, study, opts.sample, opts.lbtp,
                                    opts.tool, opts)
    print(f"{study} / {opts.sample} / {opts.lbtp}:")
    write_outputs(opts.outdir, f"{opts.sample}.{opts.lbtp}", opts.tool, tables, status)
    check_read_errors()


def run_merge(opts):
    tool = opts.tool
    if opts.file_list:
        opts.files = list(opts.files) + data_lines(opts.file_list, comment=None)
    if not opts.files and not opts.manifest:
        sys.exit("merge: no input files")
    by_table = {}
    for path in opts.files:
        table = table_of(path, tool)
        if table is None:
            sys.exit(f"{path}: not a <name>.{tool}.<table>.csv[.gz] file")
        by_table.setdefault(table, []).append(path)
    unknown = sorted(set(by_table) - set(COLLECTORS[tool][1]) - {"samples"})
    if unknown:
        sys.exit(f"tables {', '.join(unknown)} are not {tool} tables")

    os.makedirs(opts.outdir, exist_ok=True)
    print(f"{opts.name} ({tool}): merging {len(opts.files)} file(s)")

    # tables: stream rows through, so a cohort-wide merge never holds them all
    for table in COLLECTORS[tool][1]:
        paths = sorted(by_table.get(table, []))
        headers = []
        for p in paths:
            with open_text(p, "r") as fh:
                headers.append(next(csv.reader(fh), []))
        columns = table_columns_for(tool, table, headers)
        out = os.path.join(opts.outdir, f"{opts.name}.{tool}.{table}.csv.gz")
        n = 0
        with open_text(out, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
            w.writeheader()
            for p in paths:
                with open_text(p, "r") as src:
                    for row in csv.DictReader(src):
                        w.writerow(row)
                        n += 1
        print(f"  {os.path.basename(out)}: {n} row(s)")

    # samples.csv, plus a row for every manifest sample that produced nothing
    columns = ID_COLUMNS + COLLECTORS[tool][2]
    out = os.path.join(opts.outdir, f"{opts.name}.{tool}.samples.csv")
    seen, n, absent = set(), 0, 0
    with open_text(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        w.writeheader()
        for p in sorted(by_table.get("samples", [])):
            with open_text(p, "r") as src:
                for row in csv.DictReader(src):
                    seen.add((row.get("sample"), row.get("lbtp")))
                    w.writerow(row)
                    n += 1
        if opts.manifest:
            study = as_study(opts.study or opts.name)
            _, rows = manifest_rows(opts.manifest, study_code(study),
                                    opts.study_column, opts.library_types)
            for lbtp, sample in rows:
                if (sample, lbtp) not in seen:
                    w.writerow({"study": study, "sample": sample, "lbtp": lbtp,
                                "status": f"no results (no {tool} folder, or its "
                                          f"compile task failed)"})
                    absent += 1
    print(f"  {os.path.basename(out)}: {n + absent} row(s)"
          + (f", {absent} sample(s) with no results" if absent else ""))


def run_study(opts):
    study = as_study(opts.study)
    code = study_code(study)
    if opts.manifest:
        studies, samples = manifest_rows(opts.manifest, code, opts.study_column,
                                         opts.library_types)
        if not samples:
            sys.exit(f"no row of {opts.manifest} belongs to {code} "
                     f"(its studies: {', '.join(studies) or 'none'})")
    else:
        lbtps, samples = find_samples(opts.root, code, opts.library_types)
        if not samples:
            sys.exit(f"no sample folder containing '-{code}-' under {opts.root} "
                     f"(searched: {', '.join(lbtps) or 'no LBTP-* folders'})")
    print(f"{study}: {len(samples)} sample(s) in "
          f"{len({s[0] for s in samples})} library type(s)")

    def one(item):
        lbtp, sample = item
        return compile_sample(os.path.join(opts.root, lbtp, sample), study, sample,
                              lbtp, opts.tool, opts)

    merged, status_rows = {}, []
    with ThreadPoolExecutor(max_workers=max(1, opts.threads)) as pool:
        for tables, status in pool.map(one, samples):      # keeps sample order
            for t, rows in tables.items():
                merged.setdefault(t, []).extend(rows)
            status_rows += status
    write_outputs(opts.outdir, study, opts.tool, merged, status_rows)
    check_read_errors()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    modes = ap.add_subparsers(dest="mode", required=True)

    def common(p, filters=True):
        p.add_argument("--tool", required=True, choices=sorted(COLLECTORS))
        p.add_argument("--outdir", default=".")
        if filters:
            p.add_argument("--loci", help="immunotyper: file of loci (ighv, igkv, ...) "
                                          "to keep, one per line - a sample lacking "
                                          "one is reported in samples.csv")
            p.add_argument("--genes", help="immunotyper: file of genes to keep, one "
                                           "per line")

    def manifest_args(p):
        p.add_argument("--manifest", help="CSV with columns sample, lbtp[, study]")
        p.add_argument("--study-column", default="study",
                       help="manifest column giving each sample's study; where it "
                            "is absent or empty the study comes from the sample "
                            "ID (default: %(default)s)")
        p.add_argument("--library-types", nargs="*", metavar="LBTP",
                       help="only these LBTP-* folders (default: all)")

    s = modes.add_parser("sample", help="one sample's staged tool folder(s)")
    common(s)
    s.add_argument("--dir", default=".", help="folder holding the tool folder(s) "
                                              "(default: the current one)")
    s.add_argument("--study", required=True)
    s.add_argument("--sample", required=True)
    s.add_argument("--lbtp", required=True)

    m = modes.add_parser("merge", help="concatenate per-sample or per-study outputs")
    common(m, filters=False)
    m.add_argument("--name", required=True,
                   help="output name: a study (STDY-ABC) or e.g. all_cohorts")
    m.add_argument("--study", help="with --manifest: the study whose rows to "
                                   "check for (default: --name)")
    manifest_args(m)
    m.add_argument("--file-list", help="file naming the inputs, one per line - for "
                                       "more than a command line holds")
    m.add_argument("files", nargs="*", help="<name>.<tool>.<table>.csv[.gz] files")

    st = modes.add_parser("study", help="a whole study from a local results folder")
    common(st)
    st.add_argument("--root", required=True,
                    help="the results folder holding the LBTP-* folders")
    st.add_argument("--study", required=True, help="e.g. STDY-ABC (or just ABC)")
    manifest_args(st)
    st.add_argument("--threads", type=int, default=16,
                    help="samples read in parallel (default: %(default)s)")

    opts = ap.parse_args()
    if opts.mode != "merge":
        opts.loci = {x.lower() for x in names_in(opts.loci)}
        opts.genes = names_in(opts.genes)
    {"sample": run_sample, "merge": run_merge, "study": run_study}[opts.mode](opts)


if __name__ == "__main__":
    main()
