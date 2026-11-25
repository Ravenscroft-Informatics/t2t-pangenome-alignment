#!/usr/bin/env python3
import argparse, gzip, os, re, sys, subprocess, shutil
from glob import glob

try:
    import pysam
except Exception:
    pysam = None

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_line(line):
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        raise ValueError(f"Malformed line (expected ≥4 cols): {line!r}")
    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
    depth = parts[-1]
    return chrom, start, end, depth

def iter_rows_stream(path, chrom_q, qs, qe):
    # Slow fallback: scan whole file, yield overlaps
    with open_any(path) as fh:
        for line in fh:
            if not line or line[0] == "#" or line.startswith("chrom"):
                continue
            chrom, s, e, depth = parse_line(line)
            if chrom != chrom_q:           # fast reject
                continue
            if (s < qe) and (e > qs):      # half-open overlap
                yield chrom, s, e, depth

def iter_rows_tabix(path, chrom_q, qs, qe):
    # Fast path: random access via tabix
    if pysam is None:
        return None
    tbi_path = path + ".tbi"
    if not os.path.exists(tbi_path):
        return None
    tbx = pysam.TabixFile(path)  # expects bgzip + .tbi
    try:
        for line in tbx.fetch(chrom_q, qs, qe):
            # tabix does not return header lines here
            yield parse_line(line)
    except ValueError as e:
        # region not present / chrom mismatch in index
        return
    finally:
        tbx.close()

def sample_from_filename(path):
    base = os.path.basename(path)
    m = re.match(r"^(?P<sample>.+?)\.[^.]+?\.(?:per-base|regions)\.bed(?:\.gz)?$", base)
    return m.group("sample") if m else base.split(".")[0]

def open_out(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")

def expand_inputs(patterns):
    files = []
    for pat in patterns:
        expanded = glob(pat)
        if not expanded and os.path.isfile(pat):
            expanded = [pat]
        files.extend(expanded)
    return sorted(files, key=lambda p: sample_from_filename(p))

def parse_interval(interval_str, end_inclusive=False):
    s = interval_str.strip().replace(",", "")
    m = re.match(r"^([^:\s]+):(\d+)-(\d+)$", s)
    if not m:
        sys.exit(f"[ERROR] Could not parse interval '{interval_str}'. Expected 'chr:start-end'.")
    chrom, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    if end_inclusive:
        end += 1  # convert closed to half-open
    if end <= start:
        sys.exit("[ERROR] Interval end must be > start.")
    return chrom, start, end

def have_tool(name):
    return shutil.which(name) is not None

def bgzip_and_tabix(path, skip_lines=0, force=False, threads=0):
    """Create bgzip-compressed file (in place) and .tbi index if missing."""
    if not have_tool("bgzip") or not have_tool("tabix"):
        sys.stderr.write("[WARN] bgzip/tabix not found on PATH; skipping indexing.\n")
        return
    gz_path = path if path.endswith(".gz") else path + ".gz"
    if (not path.endswith(".gz")) or force:
        # bgzip in place to .gz (won't delete orig unless force=True)
        cmd = ["bgzip"]
        if threads and threads > 0:
            cmd += ["-@",
                    str(threads)]
        if force:
            cmd += ["-f"]
        cmd += ["-c", path]
        with open(gz_path, "wb") as outfh:
            subprocess.check_call(cmd, stdout=outfh)
        if force and os.path.exists(path):
            os.remove(path)
    # tabix index
    tbi = gz_path + ".tbi"
    if force or not os.path.exists(tbi):
        cmd = ["tabix", "-s", "1", "-b", "2", "-e", "3"]
        if skip_lines and skip_lines > 0:
            cmd += ["-S", str(skip_lines)]
        cmd += [gz_path]
        subprocess.check_call(cmd)
    return gz_path

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Subset mosdepth per-base/regions coverage for a target INTERVAL across many samples.\n"
            "Uses tabix for fast random access when available, falls back to streaming otherwise.\n"
            "Output: LONG format (gene, chr, start, end, sample, depth)."
        )
    )
    ap.add_argument("gene", help="Gene label to include in output (no lookup performed).")
    ap.add_argument("interval", help="Interval like 'chr11:126298670-126309536' (commas OK).")
    ap.add_argument("inputs", nargs="+", help="Input BED/BED.GZ (globs OK).")
    ap.add_argument("-o", "--out", default="coverage_subset_long.tsv.gz",
                    help="Output TSV/TSV.GZ (default: coverage_subset_long.tsv.gz)")
    ap.add_argument("--end-inclusive", action="store_true",
                    help="Treat provided end as inclusive (convert to exclusive internally).")
    ap.add_argument("--index", action="store_true",
                    help="Auto-bgzip + tabix any unindexed inputs for fast fetching.")
    ap.add_argument("--skip-lines", type=int, default=1,
                    help="Lines to skip when tabix-indexing (mosdepth often has 1 header line starting with 'chrom').")
    ap.add_argument("--threads", type=int, default=0,
                    help="Threads for bgzip during --index (passed to bgzip -@).")
    ap.add_argument("--quiet", action="store_true", help="Reduce per-file progress messages.")
    args = ap.parse_args()

    chrom_q, qs, qe = parse_interval(args.interval, end_inclusive=args.end_inclusive)
    files = expand_inputs(args.inputs)
    if not files:
        sys.exit("No input files matched.")

    # Optionally index inputs for fast random access
    if args.index:
        for i, p in enumerate(files, 1):
            gz = p if p.endswith(".gz") else p + ".gz"
            tbi = gz + ".tbi"
            if not os.path.exists(tbi):
                if not args.quiet:
                    print(f"[index] {i}/{len(files)}: {os.path.basename(p)}")
                try:
                    newp = bgzip_and_tabix(p, skip_lines=args.skip_lines, force=not p.endswith(".gz"), threads=args.threads)
                    if newp and newp != p:
                        # Replace in list if we created a .gz
                        files[files.index(p)] = newp
                except subprocess.CalledProcessError as e:
                    sys.stderr.write(f"[WARN] Failed to index {p}: {e}\n")

    n_rows_out = 0
    with open_out(args.out) as out:
        out.write("\t".join(["gene", "chr", "start", "end", "sample", "depth"]) + "\n")

        for idx, path in enumerate(files, 1):
            sample = sample_from_filename(path)
            if not args.quiet:
                print(f"[proc] {idx}/{len(files)}: {os.path.basename(path)}", file=sys.stderr)

            used_fast = False
            # Try fast tabix path first
            if pysam is not None and path.endswith(".gz") and os.path.exists(path + ".tbi"):
                for chrom, s, e, depth in iter_rows_tabix(path, chrom_q, qs, qe) or []:
                    out.write(f"{args.gene}\t{chrom}\t{s}\t{e}\t{sample}\t{depth}\n")
                    n_rows_out += 1
                used_fast = True

            if not used_fast:
                # Fallback: stream and filter (slow)
                for chrom, s, e, depth in iter_rows_stream(path, chrom_q, qs, qe):
                    out.write(f"{args.gene}\t{chrom}\t{s}\t{e}\t{sample}\t{depth}\n")
                    n_rows_out += 1

    print(f"Wrote {args.out} ({n_rows_out} rows) from {len(files)} files for {args.gene} in {chrom_q}:{qs}-{qe}")

if __name__ == "__main__":
    main()
