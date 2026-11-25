#!/Users/00104561/miniconda3/envs/mb/bin/ python
import argparse
import re
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import gzip

INDEL_SIZE_LIMIT = 40
COLORS = {"bwa": "#DE3163", "giraffe": "#44bb99"}
LABELS = {"bwa": "BWA-MEM", "giraffe": "VG-GIRAFFE"}

def mean(x):
    return float(sum(x)) / float(len(x)) if x else 0.0

def sem(x):
    return float(stats.sem(x)) if len(x) >= 2 else 0.0

def is_het(gt: str) -> bool:
    alleles = re.split(r"[\/|]", gt)
    return len(alleles) == 2 and alleles[0] != alleles[1] and alleles[0] != "." and alleles[1] != "."

def bucket_indel_len(ref: str, alts: str) -> int:
    alt0 = alts.split(",")[0]
    indel_len = len(alt0) - len(ref)
    if indel_len >= INDEL_SIZE_LIMIT:
        return INDEL_SIZE_LIMIT
    if indel_len <= -INDEL_SIZE_LIMIT:
        return -INDEL_SIZE_LIMIT
    return indel_len

def ad_to_ref_alt(ad: str):
    """Return (ref_depth, alt_depth_sum) or None."""
    if ad in (".", "", None):
        return None
    parts = ad.split(",")
    if len(parts) < 2:
        return None
    try:
        ref = 0 if parts[0] == "." else int(parts[0])
        alt = 0
        for p in parts[1:]:
            alt += 0 if p == "." else int(p)
        return ref, alt
    except ValueError:
        return None

def open_vcf(path: str):
    """Open .vcf or .vcf.gz as a text stream."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def get_fractions(vcf_path: str):
    all_counts = {"giraffe": {}, "bwa": {}}

    with open_vcf(vcf_path) as f:
        for line in f:
            if not line or line[0] == "#":
                continue

            toks = line.rstrip("\n").split("\t")
            if len(toks) < 11:
                continue  # need FORMAT + 2 samples

            ref = toks[3]
            alts = toks[4]
            fmt_keys = toks[8].split(":")

            # Assumes sample order: giraffe then bwa
            giraffe_vals = toks[9].split(":")
            bwa_vals     = toks[10].split(":")

            giraffe = {fmt_keys[i]: giraffe_vals[i] for i in range(min(len(fmt_keys), len(giraffe_vals)))}
            bwa     = {fmt_keys[i]: bwa_vals[i]     for i in range(min(len(fmt_keys), len(bwa_vals)))}

            if "GT" not in giraffe or "GT" not in bwa:
                continue
            if "AD" not in giraffe or "AD" not in bwa:
                continue
            if not (is_het(giraffe["GT"]) and is_het(bwa["GT"])):
                continue

            indel_size = bucket_indel_len(ref, alts)

            for mapper, sample_dict in (("giraffe", giraffe), ("bwa", bwa)):
                rc = ad_to_ref_alt(sample_dict.get("AD", "."))
                if rc is None:
                    continue
                ref_c, alt_c = rc
                if ref_c + alt_c == 0:
                    continue

                frac = alt_c / float(ref_c + alt_c)

                if indel_size not in all_counts[mapper]:
                    all_counts[mapper][indel_size] = [0, 0, []]
                all_counts[mapper][indel_size][0] += ref_c
                all_counts[mapper][indel_size][1] += alt_c
                all_counts[mapper][indel_size][2].append(frac)

    out = {}
    for mapper in all_counts:
        tuples = []
        for k, (ref_sum, alt_sum, fracs) in all_counts[mapper].items():
            denom = ref_sum + alt_sum
            agg = alt_sum / float(denom) if denom > 0 else 0.0
            tuples.append((k, agg, mean(fracs), sem(fracs), len(fracs)))
        out[mapper] = sorted(tuples, key=lambda x: x[0])
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf", help="2-sample merged VCF (samples: giraffe, bwa) with AD")
    ap.add_argument("out_png", help="output plot file (e.g. allele_balance.png)")
    args = ap.parse_args()

    data = get_fractions(args.vcf)

    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.size"] = 12.0
    matplotlib.rcParams["axes.grid"] = True

    fig = plt.figure(figsize=(15, 5))
    ax = plt.axes([0.07, 0.12, 0.9, 0.8])

    ax.plot([-90, 90], [0.5, 0.5], color="black", linestyle="dashed")

    for mapper in ["bwa", "giraffe"]:
        xs = [t[0] for t in data.get(mapper, [])]
        ys = [t[1] for t in data.get(mapper, [])]
        es = [t[3] for t in data.get(mapper, [])]
        ax.scatter(xs, ys, label=LABELS[mapper], color=COLORS[mapper])
        ax.plot(xs, ys, color=COLORS[mapper])
        ax.errorbar(xs, ys, yerr=es, color=COLORS[mapper], linewidth=1)

    ax.set_xlabel("Insertion or deletion length")
    ax.set_ylabel("Fraction of alternate allele")
    ax.set_xlim(-INDEL_SIZE_LIMIT - 2, INDEL_SIZE_LIMIT + 2)
    ax.set_xticks([-40, -30, -20, -10, 0, 10, 20, 30, 40])
    ax.set_xticklabels(["<-40", "-30", "-20", "-10", "0", "10", "20", "30", ">40"])
    ax.legend(loc="lower left", frameon=False)

    plt.savefig(args.out_png, dpi=200)

if __name__ == "__main__":
    main()
