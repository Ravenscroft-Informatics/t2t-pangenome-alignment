#!/usr/bin/env python3
#This script merges multiple mosdepth “regions” (or per-base) BED files into a single wide TSV, one row per genomic interval, with one column per sample and an added interval size column.
import argparse, gzip, os, re, sys
from glob import glob

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_line(line):
    # Expect: chr  start  end  gene  depth
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 5:
        raise ValueError(f"Malformed line (expected ≥5 columns): {line!r}")
    chrom, start, end, gene, depth = parts[0], parts[1], parts[2], parts[3], parts[4]
    return chrom, start, end, gene, depth

def sample_from_filename(path):
    base = os.path.basename(path)
    m = re.match(r"^(?P<sample>.+?)\.hg\d+\.(?:regions|per-base)\.bed(?:\.gz)?$", base)
    if not m:
        # Fallback: take everything before first dot
        return base.split(".")[0]
    return m.group("sample")

def read_file(path):
    with open_any(path) as fh:
        for line in fh:
            if not line or line.startswith("#") or line.startswith("chrom"):
                continue
            yield parse_line(line)

def main():
    ap = argparse.ArgumentParser(description="Aggregate mosdepth regions into a wide table.")
    ap.add_argument("inputs", nargs="+", help="Mosdepth region files (*.hg38.regions.bed[.gz])")
    ap.add_argument("-o", "--out", default="mosdepth_coverage_merged.tsv", help="Output TSV")
    args = ap.parse_args()

    files = []
    for pat in args.inputs:
        expanded = glob(pat)
        if not expanded and os.path.isfile(pat):
            expanded = [pat]
        files.extend(expanded)
    if not files:
        sys.exit("No input files matched.")
    # Sort by sample name for a stable column order
    files = sorted(files, key=lambda p: sample_from_filename(p))

    # Read the first file as the reference intervals
    ref_path = files[0]
    ref_rows = []
    for row in read_file(ref_path):
        ref_rows.append(row[:4])  # (chr, start, end, gene)

    # Prepare per-sample depth columns; verify intervals match
    sample_names = []
    depths_by_sample = []  # list of lists aligned to ref_rows
    for idx, path in enumerate(files):
        sname = sample_from_filename(path)
        sample_names.append(sname)
        depths = []
        it = read_file(path)
        for i, row in enumerate(it):
            key = row[:4]
            if i >= len(ref_rows):
                sys.exit(f"[ERROR] {os.path.basename(path)} has MORE rows than reference {os.path.basename(ref_path)} (first extra at line {i+1}: {key})")
            if key != ref_rows[i]:
                rkey = ref_rows[i]
                sys.exit(
                    f"[ERROR] Interval mismatch at line {i+1}.\n"
                    f"  Reference ({os.path.basename(ref_path)}): {rkey}\n"
                    f"  This file  ({os.path.basename(path)}): {key}"
                )
            depths.append(row[4])
        if len(depths) != len(ref_rows):
            sys.exit(f"[ERROR] {os.path.basename(path)} has FEWER rows ({len(depths)}) than reference ({len(ref_rows)}).")
        depths_by_sample.append(depths)

    # Write output
    with open(args.out, "w") as out:
        header = ["chr", "start", "end", "gene"] + sample_names
        out.write("\t".join(header) + "\n")
        nrows = len(ref_rows)
        nsamples = len(sample_names)
        for i in range(nrows):
            row_key = ref_rows[i]
            values = [depths_by_sample[s][i] for s in range(nsamples)]
            out.write("\t".join(row_key + tuple(values)) + "\n")

    print(f"Wrote {args.out} with {len(ref_rows)} intervals and {len(sample_names)} samples.")

if __name__ == "__main__":
    main()
