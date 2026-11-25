#!/usr/bin/env python3
"""
Classify vg giraffe GAF alignments by whether they stay on the reference path
(e.g. GRCh38) or use sequence/edges that exist only in non-reference haplotypes,
and report alignment QC (perfect reads, MAPQ 0, proper pairing, mismatch rate).

Inputs
  --ref-gaf  : the reference path exported from the SAME .gbz used for giraffe:
               vg paths -x graph.gbz -S GRCh38 -A > grch38.paths.gaf
  --reads    : giraffe GAF (plain or .gz), node-ID path field (default giraffe output,
               i.e. NOT run with --named-coordinates)

Per-read categories (primary alignments only; giraffe GAF is primary-only by default)
  unmapped      path == '*'
  low_mapq      mapped, MAPQ < --min-mapq
  ref_only      every node is on the reference path AND every node-to-node step
                follows the reference path (read is fully explained by GRCh38)
  nonref_edge   every node is on the reference path but at least one step is
                not a reference edge (deletion / rearrangement vs GRCh38)
  nonref_node   at least one node is absent from the reference path
                (insertion / SNV allele / divergent sequence from other haplotypes)

Per-category metrics (table)
  perfect_reads       whole read aligned (no soft-clip), no mismatches/indels vs its
                      graph path: qstart==0, qend==qlen, matches==block==qlen
  mapq0_reads         MAPQ == 0
  proper_pair_reads   giraffe pd:b:1 tag
  mismatch_rate       mismatches / (matches + mismatches), from the cs tag
  error_rate          (mismatches + inserted + deleted bases) / (matches + mismatches
                      + inserted bases) -- same definition as samtools stats
                      "error rate" (NM / bases mapped (cigar))
  NB: mismatches/errors are against the read's GRAPH PATH, not GRCh38.

Summary metrics over ALL mapped reads (classified + low_mapq) and pair counts
(in reads, no MAPQ filter) are written as '## key<TAB>value' lines.

Usage
  python3 gaf_ref_vs_alt.py --ref-gaf grch38.paths.gaf --reads sample.gaf.gz \
      [--min-mapq 0] [--nonref-names nonref_reads.txt.gz] > sample.ref_vs_alt.tsv

Faster rerun (--known-nonref)
  Pass the nonref_reads.txt.gz written by an earlier run on the SAME GAF (with
  --min-mapq <= the current value). Reads listed there take their category from the
  file; every other mapped read is ref_only, so path classification is only done for
  nonref_node reads (to count their non-ref node visits). The file is streamed in step
  with the GAF (both are in GAF order), so it costs no extra memory. The run aborts if
  the file does not match the GAF. Reading the GAF from a separate decompressor also
  helps: pigz -dc sample.gaf.gz | python3 gaf_ref_vs_alt.py --reads - ...
"""
import argparse
import gzip
import re
import sys

import numpy as np

STEP_RE = re.compile(r"([<>])(\d+)")
CS_RE = re.compile(r"(:)(\d+)|(\*)[A-Za-z][A-Za-z]|([+\-=])([A-Za-z]+)")


def xopen(path, mode="rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def load_reference(ref_gaf):
    """Return (rank, orient).

    rank[node]   : global position of the node along the concatenated reference
                   paths (gap of 2 between contigs so contig ends never look adjacent);
                   -1 = not on reference; -2 = on reference more than once (ambiguous).
    orient[node] : +1 if the reference traverses the node forward, -1 if reverse.
    """
    contigs = []
    max_id = 0
    with xopen(ref_gaf) as fh:
        for line in fh:
            f = line.split("\t", 6)
            if len(f) < 6 or f[5] == "*":
                continue
            steps = STEP_RE.findall(f[5])
            ids = np.fromiter((int(n) for _, n in steps), dtype=np.int64, count=len(steps))
            ori = np.fromiter((1 if o == ">" else -1 for o, _ in steps), dtype=np.int8, count=len(steps))
            contigs.append((ids, ori))
            if len(ids):
                max_id = max(max_id, int(ids.max()))
            print(f"[ref] {f[0]}: {len(ids):,} nodes", file=sys.stderr)

    rank = np.full(max_id + 1, -1, dtype=np.int64)
    orient = np.zeros(max_id + 1, dtype=np.int8)
    pos = 0
    for ids, ori in contigs:
        r = np.arange(pos, pos + len(ids), dtype=np.int64)
        seen = rank[ids] != -1
        rank[ids[seen]] = -2  # also on an earlier contig -> ambiguous
        fresh = ids[~seen]
        rank[fresh] = r[~seen]
        orient[fresh] = ori[~seen]
        uniq, counts = np.unique(ids, return_counts=True)
        rank[uniq[counts > 1]] = -2  # repeated within this contig
        pos += len(ids) + 2
    n_ref = int((rank != -1).sum())
    n_amb = int((rank == -2).sum())
    print(f"[ref] {n_ref:,} reference nodes ({n_amb:,} visited >1x, treated as ref, edge check relaxed)",
          file=sys.stderr)
    return rank, orient


def classify(path, rank, orient, nmax):
    """Return (category, n_nodes, n_nonref_nodes)."""
    steps = STEP_RE.findall(path)
    n_nonref = 0
    nonref_edge = False
    prev_rank = prev_strand = None
    for o, n in steps:
        nid = int(n)
        r = rank[nid] if nid < nmax else -1
        if r == -1:
            n_nonref += 1
            prev_rank = None
            continue
        if r == -2:  # ambiguous repeat node: on reference, can't check adjacency
            prev_rank = None
            continue
        strand = (1 if o == ">" else -1) * orient[nid]
        if prev_rank is not None and not nonref_edge:
            if strand != prev_strand or r - prev_rank != strand:
                nonref_edge = True
        prev_rank, prev_strand = r, strand
    if n_nonref:
        return "nonref_node", len(steps), n_nonref
    if nonref_edge:
        return "nonref_edge", len(steps), 0
    return "ref_only", len(steps), 0


def cs_counts(cs):
    """Return (matches, mismatches, inserted_bases, deleted_bases) from a cs string."""
    m = mm = ins = dele = 0
    for colon, num, star, op, seq in CS_RE.findall(cs):
        if colon:
            m += int(num)
        elif star:
            mm += 1
        elif op == "+":
            ins += len(seq)
        elif op == "-":
            dele += len(seq)
        elif op == "=":
            m += len(seq)
    return m, mm, ins, dele


def get_tag(tags, key):
    """Value of TAG:TYPE (e.g. 'cs:Z:') from the tab-separated tag string, or None."""
    i = tags.find(key)
    while i > 0 and tags[i - 1] != "\t":
        i = tags.find(key, i + 1)
    if i < 0:
        return None
    j = tags.find("\t", i)
    return tags[i + len(key): j if j >= 0 else None]


class Stats:
    __slots__ = ("reads", "perfect", "mapq0", "proper", "m", "mm", "ins", "dele")

    def __init__(self):
        for s in self.__slots__:
            setattr(self, s, 0)

    def add(self, other):
        for s in self.__slots__:
            setattr(self, s, getattr(self, s) + getattr(other, s))
        return self


def frac(x, y):
    return x / y if y else float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref-gaf", required=True)
    ap.add_argument("--reads", required=True, help="giraffe GAF (.gaf or .gaf.gz, '-' for stdin)")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--nonref-names", help="optional output (.gz ok): name, category, MAPQ, perfect(0/1) "
                                           "of non-ref reads")
    ap.add_argument("--known-nonref", help="nonref_reads.txt(.gz) from an earlier run on the same GAF; "
                                           "skips classifying reads (see 'Faster rerun')")
    a = ap.parse_args()

    rank, orient = load_reference(a.ref_gaf)
    nmax = len(rank)

    cats = ["unmapped", "low_mapq", "ref_only", "nonref_edge", "nonref_node"]
    classified = ("ref_only", "nonref_edge", "nonref_node")
    st = {c: Stats() for c in cats}
    nodes_total = nodes_nonref = 0
    names_out = xopen(a.nonref_names, "wt") if a.nonref_names else None

    # pairing (counted in reads)
    pend = {}
    pr = dict.fromkeys(["paired_reads", "unpaired_reads", "both_mates_mapped", "mate_unmapped"], 0)

    def do_pair(x_mapped, y_mapped):
        pr["paired_reads"] += 2
        if x_mapped and y_mapped:
            pr["both_mates_mapped"] += 2
        elif x_mapped or y_mapped:
            pr["mate_unmapped"] += 1

    # --known-nonref: stream (name, category) in GAF order
    known = xopen(a.known_nonref) if a.known_nonref else None
    known_name = known_cat = None
    n_known_used = n_known_disagree = 0

    def next_known():
        nonlocal known_name, known_cat
        for kl in known:
            k = kl.rstrip("\n").split("\t", 2)
            if len(k) >= 2:
                known_name, known_cat = k[0], k[1]
                return
        known_name = known_cat = None

    if known:
        next_known()

    with xopen(a.reads) as fh:
        for i, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t", 12)
            tags = f[12] if len(f) > 12 else ""
            name = f[0]
            paired = name[-2:] in ("/1", "/2") or "fn:Z:" in tags or "fp:Z:" in tags

            if f[5] == "*":
                st["unmapped"].reads += 1
                info = False
            else:
                mapq = int(f[11])
                if known is None:
                    cat, nn, nnr = classify(f[5], rank, orient, nmax)
                elif name == known_name:
                    n_known_used += 1
                    cat = known_cat
                    if cat == "nonref_node":
                        c2, nn, nnr = classify(f[5], rank, orient, nmax)
                        n_known_disagree += c2 != cat
                    else:
                        nn, nnr = f[5].count(">") + f[5].count("<"), 0
                    next_known()
                else:
                    cat = "ref_only"
                    nn, nnr = f[5].count(">") + f[5].count("<"), 0
                info = True
                if mapq < a.min_mapq:
                    cat = "low_mapq"
                else:
                    nodes_total += nn
                    nodes_nonref += nnr
                s = st[cat]
                s.reads += 1
                if mapq == 0:
                    s.mapq0 += 1
                if "pd:b:1" in tags:
                    s.proper += 1
                qlen = f[1]
                is_perfect = f[2] == "0" and f[3] == qlen and f[9] == qlen and f[10] == qlen
                if is_perfect:
                    s.perfect += 1
                    s.m += int(qlen)
                else:
                    cs = get_tag(tags, "cs:Z:")
                    if cs is not None:
                        m, mm, ins, dele = cs_counts(cs)
                    else:  # fallback: block - matches lumps mismatches + gaps
                        m, mm, ins, dele = int(f[9]), int(f[10]) - int(f[9]), 0, 0
                    s.m += m
                    s.mm += mm
                    s.ins += ins
                    s.dele += dele
                if names_out and cat in ("nonref_edge", "nonref_node"):
                    names_out.write(f"{name}\t{cat}\t{mapq}\t{int(is_perfect)}\n")

            if paired:
                key = name[:-2] if name[-2:] in ("/1", "/2") else name
                mate = pend.pop(key, None)
                if mate is None:
                    pend[key] = info
                else:
                    do_pair(mate, info)
            else:
                pr["unpaired_reads"] += 1

            if i % 10_000_000 == 0:
                print(f"[reads] {i:,} processed ({len(pend):,} mates pending)", file=sys.stderr)
    if names_out:
        names_out.close()
    if known:
        leftover = known_name is not None
        known.close()
        print(f"[known-nonref] {n_known_used:,} reads taken from {a.known_nonref}; "
              f"{n_known_disagree:,} nonref_node disagreements", file=sys.stderr)
        if leftover or n_known_disagree:
            sys.exit(f"ERROR: {a.known_nonref} does not match {a.reads} "
                     f"({'unconsumed entries remain' if leftover else 'classification disagreements'}); "
                     f"rerun without --known-nonref")
    orphans = len(pend)
    pr["unpaired_reads"] += orphans  # mate never seen in the file

    # ------------------------------------------------------------------ output
    rows = [(c, st[c]) for c in cats]
    rows.append(("any_nonref", Stats().add(st["nonref_edge"]).add(st["nonref_node"])))
    allc = Stats()
    for c in classified:
        allc.add(st[c])
    rows.append(("all_classified", allc))
    total = sum(st[c].reads for c in cats)
    mapped_all = Stats().add(allc).add(st["low_mapq"])

    print("category\treads\tfrac_of_all\tfrac_of_classified\tperfect_reads\tfrac_perfect_in_category"
          "\tmapq0_reads\tfrac_mapq0_in_category\tproper_pair_reads\tfrac_proper_pair_in_category"
          "\tmismatch_rate\terror_rate")
    for c, s in rows:
        fm = frac(s.reads, allc.reads) if c not in ("unmapped", "low_mapq") else float("nan")
        if c == "unmapped":
            print(f"{c}\t{s.reads}\t{frac(s.reads, total):.6f}\tnan\t0\tnan\t0\tnan\t0\tnan\tnan\tnan")
            continue
        mrate = frac(s.mm, s.m + s.mm)
        erate = frac(s.mm + s.ins + s.dele, s.m + s.mm + s.ins)
        print(f"{c}\t{s.reads}\t{frac(s.reads, total):.6f}\t{fm:.6f}\t{s.perfect}\t{frac(s.perfect, s.reads):.6f}"
              f"\t{s.mapq0}\t{frac(s.mapq0, s.reads):.6f}\t{s.proper}\t{frac(s.proper, s.reads):.6f}"
              f"\t{mrate:.6g}\t{erate:.6g}")

    p_all, p_nonref = allc.perfect, rows[5][1].perfect
    print(f"# perfect reads as fraction of ALL reads (incl. unmapped/low-MAPQ): "
          f"{p_all}/{total} ({frac(p_all, total):.6f})")
    print(f"# share of perfect reads that are non-reference: {p_nonref}/{p_all} "
          f"({frac(p_nonref, p_all):.6f})")
    print(f"# node visits on non-ref nodes: {nodes_nonref}/{nodes_total} "
          f"({frac(nodes_nonref, nodes_total):.6f})")

    # all mapped reads (classified + low_mapq): key<TAB>value
    m = mapped_all
    kv = [
        ("mapped_reads", m.reads),
        ("mapped_mapq0_reads", m.mapq0),
        ("mapped_frac_mapq0", f"{frac(m.mapq0, m.reads):.6f}"),
        ("mapped_proper_pair_reads", m.proper),
        ("mapped_frac_proper_pair", f"{frac(m.proper, m.reads):.6f}"),
        ("mapped_matched_bases", m.m),
        ("mapped_mismatched_bases", m.mm),
        ("mapped_inserted_bases", m.ins),
        ("mapped_deleted_bases", m.dele),
        ("mapped_mismatch_rate", f"{frac(m.mm, m.m + m.mm):.6g}"),
        ("mapped_error_rate", f"{frac(m.mm + m.ins + m.dele, m.m + m.mm + m.ins):.6g}"),
    ] + [(k, v) for k, v in pr.items()] + [
        ("orphan_reads", orphans),
    ]
    for k, v in kv:
        print(f"## {k}\t{v}")


if __name__ == "__main__":
    main()