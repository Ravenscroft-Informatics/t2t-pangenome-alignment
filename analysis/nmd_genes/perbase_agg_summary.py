#!/usr/bin/env python3
#This script computes per-(gene, sample) summary statistics from a pre-sorted TSV of per-base (or per-interval) depths, and writes one output row per (gene, sample).
import argparse, csv, math

def idxs_for_median(n:int):
    # 0-based indices of the median position(s)
    if n % 2 == 1:
        i = n // 2
        return (i, i)
    else:
        j = n // 2
        return (j-1, j)

def main():
    ap = argparse.ArgumentParser(
        description="Exact min/max/range/mean/median/SD per (gene,sample) from a file pre-sorted by gene,sample,depth."
    )
    ap.add_argument("-i","--input", required=True, help="Sorted TSV (gene,chr,start,end,sample,depth)")
    ap.add_argument("-o","--output", required=True, help="Output TSV")
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--depth-col", default="depth")
    ap.add_argument("--stdev", choices=["sample","population"], default="sample")
    args = ap.parse_args()

    # ---------- PASS 1: counts + Welford + min/max ----------
    counts = {}  # (gene,sample) -> n
    means  = {}  # -> mean
    M2s    = {}  # -> sum of squares of diffs from mean
    mins   = {}  # -> min depth
    maxs   = {}  # -> max depth

    with open(args.input, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        try:
            gi = header.index(args.gene_col)
            si = header.index(args.sample_col)
            di = header.index(args.depth_col)
        except ValueError as e:
            raise SystemExit(f"Missing required column: {e}. Found: {header}")

        for row in reader:
            try:
                d = float(row[di])
            except Exception:
                continue
            key = (row[gi], row[si])

            n = counts.get(key, 0) + 1
            counts[key] = n

            # Welford update
            mean = means.get(key, 0.0)
            delta = d - mean
            mean += delta / n
            means[key] = mean
            delta2 = d - mean
            M2s[key] = M2s.get(key, 0.0) + delta * delta2

            # min/max (file sorted by depth ensures these converge quickly, but we don’t rely on that)
            mins[key] = d if key not in mins else min(mins[key], d)
            maxs[key] = d if key not in maxs else max(maxs[key], d)

    # Precompute median target indices for each key
    med_idx = {k: idxs_for_median(n) for k, n in counts.items()}

    # ---------- PASS 2: pick medians from the depth-sorted stream ----------
    # We track per-key row index within the (gene,sample) group
    seen = {}  # (gene,sample) -> index within group (0-based)
    med_left  = {}  # -> left median value
    med_right = {}  # -> right median value

    with open(args.input, "r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        _ = next(reader)  # skip header

        current_key = None
        for row in reader:
            try:
                d = float(row[di])
            except Exception:
                continue
            key = (row[gi], row[si])

            # group boundary? reset counter
            if key != current_key:
                current_key = key
                seen[key] = 0
            else:
                seen[key] += 1

            i = seen[key]
            a, b = med_idx[key]
            if i == a:
                med_left[key] = d
            if i == b:
                med_right[key] = d

    # ---------- Write results ----------
    with open(args.output, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene","sample","min_depth","max_depth","range_depth",
                    "mean_depth","median_depth","sd_depth"])

        for key in counts:
            n = counts[key]
            mean = means[key]
            var = (M2s[key] / (n-1)) if (args.stdev=="sample" and n>1) else ((M2s[key] / n) if n>0 else 0.0)
            sd = math.sqrt(var) if n>1 else 0.0
            mn = mins[key]; mx = maxs[key]
            rng = mx - mn
            med = (med_left[key] + med_right[key]) / 2.0
            gene, sample = key
            w.writerow([gene, sample, f"{mn:g}", f"{mx:g}", f"{rng:g}",
                        f"{mean:g}", f"{med:g}", f"{sd:g}"])

if __name__ == "__main__":
    main()
