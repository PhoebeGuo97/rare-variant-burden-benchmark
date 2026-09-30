#!/usr/bin/env python3
"""
Day 1, step 1b: measure segmental duplication overlap per gene.

This replaces the earlier hand-assertion that certain genes "have bad
mappability". Overlap with the UCSC hg38 segmental duplication track is an
objective, citable annotation. It is computed for every gene in the panel, so no
gene is singled out by assumption.

This is an a priori flag only. The empirical check comes at the slicing stage
(call rate, HWE, excess heterozygosity in the actual 1000G callset), and where
the two disagree the empirical one wins.

Inputs
  config/panel_genes_coords.tsv
  ref/genomicSuperDups.txt.gz      UCSC hg38 segmental duplications

Output
  config/gene_mappability.tsv      symbol, segdup_bp, segdup_frac, flag

Usage
  python3 scripts/01b_mappability.py
"""

import argparse
import gzip
from collections import defaultdict

# Fraction of the gene body overlapping a segmental duplication above which the
# gene is flagged. 0.5 is a convention, not a derived threshold. The continuous
# segdup_frac column is what the analysis should actually use.
FLAG_THRESHOLD = 0.5


def load_segdups(path):
    """UCSC genomicSuperDups: bin, chrom, chromStart, chromEnd, ... (0-based BED)."""
    by_chrom = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            by_chrom[f[1]].append((int(f[2]), int(f[3])))
    # merge so overlapping duplication records are not double counted
    merged = {}
    for c, ivs in by_chrom.items():
        out, cur = [], None
        for s, e in sorted(ivs):
            if cur and s <= cur[1]:
                cur = (cur[0], max(cur[1], e))
            else:
                if cur:
                    out.append(cur)
                cur = (s, e)
        if cur:
            out.append(cur)
        merged[c] = out
    return merged


def overlap_bp(segdups, chrom, start, end):
    """start/end are 1-based inclusive from the GFF; convert to 0-based half-open."""
    s0, e0 = start - 1, end
    total = 0
    for s, e in segdups.get("chr" + chrom, []):
        if e <= s0:
            continue
        if s >= e0:
            break
        total += min(e, e0) - max(s, s0)
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", default="config/panel_genes_coords.tsv")
    ap.add_argument("--segdups", default="ref/genomicSuperDups.txt.gz")
    ap.add_argument("--out", default="config/gene_mappability.tsv")
    args = ap.parse_args()

    segdups = load_segdups(args.segdups)
    print("segdup records loaded for %d chromosomes" % len(segdups))

    rows = []
    with open(args.coords) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {k: i for i, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            sym = f[idx["symbol"]]
            chrom = f[idx["chrom"]]
            start, end = int(f[idx["start"]]), int(f[idx["end"]])
            glen = end - start + 1
            bp = overlap_bp(segdups, chrom, start, end)
            frac = bp / glen
            rows.append((sym, f[idx["list"]], f[idx["cluster"]], glen, bp, frac))

    rows.sort(key=lambda r: -r[5])

    with open(args.out, "w") as out:
        out.write("symbol\tlist\tcluster\tgene_bp\tsegdup_bp\tsegdup_frac\tsegdup_flag\n")
        for sym, lst, cl, glen, bp, frac in rows:
            flag = "high" if frac >= FLAG_THRESHOLD else "low"
            out.write("%s\t%s\t%s\t%d\t%d\t%.3f\t%s\n"
                      % (sym, lst, cl, glen, bp, frac, flag))

    n_high = sum(1 for r in rows if r[5] >= FLAG_THRESHOLD)
    n_any = sum(1 for r in rows if r[5] > 0)

    print()
    print("genes with any segdup overlap : %d / %d" % (n_any, len(rows)))
    print("genes flagged high (>= %.0f%%)  : %d" % (FLAG_THRESHOLD * 100, n_high))
    print()
    print("--- top 12 by segdup fraction ---")
    print("%-12s %-16s %-12s %8s  %s" % ("symbol", "list", "cluster", "frac", "flag"))
    for sym, lst, cl, glen, bp, frac in rows[:12]:
        print("%-12s %-16s %-12s %8.3f  %s"
              % (sym, lst, cl, frac, "high" if frac >= FLAG_THRESHOLD else "low"))
    print()
    print("wrote %s" % args.out)
    print()
    print("Read this against expectation. If the 17q21.31 genes do NOT come out")
    print("high, the assumption that they are artifact prone was wrong and the")
    print("stratification should follow this table, not the prior.")


if __name__ == "__main__":
    main()
