#!/usr/bin/env python3
import argparse, gzip, os, re, sys, subprocess, shutil
from glob import glob
from collections import defaultdict

try:
    import pysam
except Exception:
    pysam = None


# -------------------- IO helpers --------------------

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def open_out(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")


def have_tool(name):
    return shutil.which(name) is not None


def bgzip_and_tabix(path, skip_lines=0, force=False, threads=0):
    """Create bgzip-compressed file (in place) and .tbi index if missing."""
    if not have_tool("bgzip") or not have_tool("tabix"):
        sys.stderr.write("[WARN] bgzip/tabix not found on PATH; skipping indexing.\n")
        return
    gz_path = path if path.endswith(".gz") else path + ".gz"
    if (not path.endswith(".gz")) or force:
        cmd = ["bgzip"]
        if threads and threads > 0:
            cmd += ["-@", str(threads)]
        if force:
            cmd += ["-f"]
        cmd += ["-c", path]
        with open(gz_path, "wb") as outfh:
            subprocess.check_call(cmd, stdout=outfh)
        if force and os.path.exists(path):
            os.remove(path)
    tbi = gz_path + ".tbi"
    if force or not os.path.exists(tbi):
        cmd = ["tabix", "-s", "1", "-b", "2", "-e", "3"]
        if skip_lines and skip_lines > 0:
            cmd += ["-S", str(skip_lines)]
        cmd += [gz_path]
        subprocess.check_call(cmd)
    return gz_path


# -------------------- parsing helpers --------------------

def parse_mosdepth_line(line):
    # Expect: chr  start  end  ... depth (last col)
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        raise ValueError(f"Malformed line (expected ≥4 cols): {line!r}")
    chrom, start, end = parts[0], int(parts[1]), int(parts[2])
    depth = parts[-1]
    return chrom, start, end, depth


def sample_from_filename(path):
    base = os.path.basename(path)
    m = re.match(r"^(?P<sample>.+?)\.[^.]+?\.(?:per-base|regions)\.bed(?:\.gz)?$", base)
    return m.group("sample") if m else base.split(".")[0]


def expand_inputs(patterns):
    files = []
    for pat in patterns:
        expanded = glob(pat)
        if not expanded and os.path.isfile(pat):
            expanded = [pat]
        files.extend(expanded)
    return sorted(files, key=lambda p: sample_from_filename(p))


def load_bed(bed_path, assume_bed_end_inclusive=False):
    """
    Load a BED/BED.GZ file with >=3 cols. Column 4 (if present) is used as the gene/region label.
    Returns a list of dicts: {'chrom','start','end','label'}.
    BED is standard: 0-based, end-exclusive. If `assume_bed_end_inclusive` is True,
    convert to half-open by adding +1 to end.
    """
    out = []
    with open_any(bed_path) as fh:
        for ln in fh:
            if not ln.strip():
                continue
            if ln.startswith("#") or ln.startswith("track") or ln.startswith("browser"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                # likely a header like "chrom start end ..."
                continue
            if assume_bed_end_inclusive:
                end += 1
            if end <= start:
                continue
            label = parts[3] if len(parts) >= 4 and parts[3] else f"{chrom}:{start}-{end}"
            out.append({"chrom": chrom, "start": start, "end": end, "label": label})
    if not out:
        sys.exit(f"[ERROR] No usable intervals found in BED: {bed_path}")
    return out


# -------------------- iterators over mosdepth tables --------------------

def iter_rows_stream(path, chrom_q, qs, qe):
    """Slow fallback: scan whole file, yield overlaps with [qs,qe) on chrom_q."""
    with open_any(path) as fh:
        for line in fh:
            if not line or line[0] == "#" or line.startswith("chrom"):
                continue
            chrom, s, e, depth = parse_mosdepth_line(line)
            if chrom != chrom_q:           # fast reject
                continue
            if (s < qe) and (e > qs):      # half-open overlap
                yield chrom, s, e, depth


def iter_rows_tabix(path, chrom_q, qs, qe):
    """Fast path: random access via tabix (requires bgzip + .tbi)."""
    if pysam is None:
        return None
    tbi_path = path + ".tbi"
    if not os.path.exists(tbi_path):
        return None
    tbx = pysam.TabixFile(path)
    try:
        for line in tbx.fetch(chrom_q, qs, qe):
            yield parse_mosdepth_line(line)
    except ValueError:
        # region not present / chrom mismatch in index
        return
    finally:
        tbx.close()


# -------------------- main --------------------

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Subset mosdepth per-base/regions coverage for ALL intervals in a BED file across many samples.\n"
            "Uses tabix for fast random access when available, falls back to streaming otherwise.\n"
            "Output: LONG format (gene, chr, start, end, sample, depth)."
        )
    )
    ap.add_argument("bed",
                    help="BED/BED.GZ with >=3 columns (chrom, start, end). "
                         "If column 4 exists, it's used as the gene/region label.")
    ap.add_argument("inputs", nargs="+",
                    help="Input mosdepth per-base/regions BED/BED.GZ (globs OK).")
    ap.add_argument("-o", "--out", default="coverage_subset_long.tsv.gz",
                    help="Output TSV/TSV.GZ (default: coverage_subset_long.tsv.gz)")
    ap.add_argument("--bed-end-inclusive", action="store_true",
                    help="Treat BED 'end' as inclusive (convert to exclusive internally).")
    ap.add_argument("--index", action="store_true",
                    help="Auto-bgzip + tabix any unindexed inputs for fast fetching.")
    ap.add_argument("--skip-lines", type=int, default=1,
                    help="Lines to skip when tabix-indexing (mosdepth often has 1 header line starting with 'chrom').")
    ap.add_argument("--threads", type=int, default=0,
                    help="Threads for bgzip during --index (passed to bgzip -@).")
    ap.add_argument("--quiet", action="store_true",
                    help="Reduce per-file progress messages.")
    args = ap.parse_args()

    intervals = load_bed(args.bed, assume_bed_end_inclusive=args.bed_end_inclusive)

    # Group intervals by chromosome for a tiny speed-up with tabix
    intervals_by_chr = defaultdict(list)
    for iv in intervals:
        intervals_by_chr[iv["chrom"]].append(iv)

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
                    newp = bgzip_and_tabix(
                        p, skip_lines=args.skip_lines,
                        force=not p.endswith(".gz"),
                        threads=args.threads
                    )
                    if newp and newp != p:
                        files[files.index(p)] = newp
                except subprocess.CalledProcessError as e:
                    sys.stderr.write(f"[WARN] Failed to index {p}: {e}\n")

    n_rows_out = 0
    with open_out(args.out) as out:
        out.write("\t".join(["gene", "chr", "start", "end", "sample", "depth"]) + "\n")

        for f_idx, path in enumerate(files, 1):
            sample = sample_from_filename(path)
            if not args.quiet:
                print(f"[proc] {f_idx}/{len(files)}: {os.path.basename(path)}", file=sys.stderr)

            fast_ok = (pysam is not None and path.endswith(".gz") and os.path.exists(path + ".tbi"))

            if fast_ok:
                # Tabix path: fetch per interval
                for chrom_q, ivs in intervals_by_chr.items():
                    for iv in ivs:
                        it = iter_rows_tabix(path, chrom_q, iv["start"], iv["end"])
                        if it is None:
                            # fallback for this file if somehow tabix couldn't fetch
                            for chrom, s, e, depth in iter_rows_stream(path, chrom_q, iv["start"], iv["end"]):
                                out.write(f"{iv['label']}\t{chrom}\t{s}\t{e}\t{sample}\t{depth}\n")
                                n_rows_out += 1
                            continue
                        for chrom, s, e, depth in it:
                            # Half-open overlap is already enforced by tabix fetch range
                            out.write(f"{iv['label']}\t{chrom}\t{s}\t{e}\t{sample}\t{depth}\n")
                            n_rows_out += 1
            else:
                # Streaming path: iterate whole file once and test overlap with *any* intervals on that chrom
                # To keep memory modest, build an interval list per chrom; we already have that.
                with open_any(path) as fh:
                    for line in fh:
                        if not line or line[0] == "#" or line.startswith("chrom"):
                            continue
                        chrom, s, e, depth = parse_mosdepth_line(line)
                        if chrom not in intervals_by_chr:
                            continue
                        for iv in intervals_by_chr[chrom]:
                            if (s < iv["end"]) and (e > iv["start"]):
                                out.write(f"{iv['label']}\t{chrom}\t{s}\t{e}\t{sample}\t{depth}\n")
                                n_rows_out += 1

    print(f"Wrote {args.out} ({n_rows_out} rows) from {len(files)} files and {len(intervals)} BED intervals")

if __name__ == "__main__":
    main()
