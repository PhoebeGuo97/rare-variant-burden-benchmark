#!/usr/bin/env python3
"""
Day 2, step 7b: turn bcftools csq output into REGENIE mask files.

Inputs
  results/csq_table.tsv.gz          ID, CHROM, POS, BCSQ   (from 07_make_masks.sh)
  results/variant_qc_pooled.tsv.gz  AC, AN                 (from 04_variant_qc.sh)
  config/panel_genes_coords.tsv     the 56 panel genes

Outputs
  config/regenie_annot.tsv     varID  gene  annotation
  config/regenie_setlist.tsv   gene  chrom  pos  comma-separated varIDs
  config/regenie_masks.def     mask_name  annotation[,annotation...]
  results/mask_counts.tsv      testable variants per gene x class x AAF bin

Two decisions that change the answer, both stated rather than buried.

1. ONLY protein_coding transcripts count. The BCSQ biotype field also carries
   NMD, retained_intron and similar. A stop_gained on an NMD-target transcript
   is not a loss-of-function event for the gene, and counting it would inflate
   the strictest mask with exactly the variants it is supposed to be pure of.

2. Each variant takes its MOST SEVERE consequence across the gene's
   protein_coding transcripts. The alternative is to use one canonical or MANE
   transcript. Most-severe is the more inclusive convention and is what most
   published burden pipelines use, but it does raise LoF counts relative to a
   canonical-only definition, so it is a real choice and not a default.
"""

import argparse
import gzip
import sys
from collections import defaultdict

LOF = {"stop_gained", "frameshift", "splice_acceptor", "splice_donor",
       "start_lost", "stop_lost"}
MISSENSE = {"missense", "inframe_deletion", "inframe_insertion",
            "inframe_altering"}
SYNONYMOUS = {"synonymous", "stop_retained", "start_retained"}

# severity order, most severe first
CLASS_RANK = ["LoF", "missense", "synonymous"]

MASKS = [
    ("M1_lof", ["LoF"]),
    ("M2_lof_mis", ["LoF", "missense"]),
    ("M4_all", ["LoF", "missense", "synonymous"]),
]

AAF_BINS = [0.001, 0.005, 0.01]


def classify(term):
    if term in LOF:
        return "LoF"
    if term in MISSENSE:
        return "missense"
    if term in SYNONYMOUS:
        return "synonymous"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csq", default="results/csq_table.tsv.gz")
    ap.add_argument("--freq", default="results/variant_qc_pooled.tsv.gz")
    ap.add_argument("--coords", default="config/panel_genes_coords.tsv")
    ap.add_argument("--outdir", default="config")
    args = ap.parse_args()

    # panel genes and their set-list anchor position
    panel, anchor = set(), {}
    with open(args.coords) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        i = {k: n for n, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(h):
                continue
            panel.add(f[i["symbol"]])
            anchor[f[i["symbol"]]] = (f[i["chrom"]], int(f[i["start"]]))
    print("panel genes: %d" % len(panel))

    # allele frequency, keyed by the same chr:pos:ref:alt the IDs use
    maf = {}
    with gzip.open(args.freq, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(h):
                continue
            try:
                ac, an = float(f[4]), float(f[5])
            except ValueError:
                continue
            if an <= 0:
                continue
            af = ac / an
            maf["%s:%s:%s:%s" % (f[0], f[1], f[2], f[3])] = min(af, 1 - af)
    print("variants with a frequency: %d" % len(maf))

    # variant -> gene -> best class
    assign = defaultdict(dict)
    n_rows = n_with_csq = n_nonpanel = 0
    seen_biotypes = defaultdict(int)
    dropped_biotype = 0

    with gzip.open(args.csq, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            n_rows += 1
            vid, bcsq = f[0], f[3]
            if bcsq in (".", ""):
                continue
            n_with_csq += 1
            for entry in bcsq.split(","):
                entry = entry.lstrip("@*")
                p = entry.split("|")
                if len(p) < 4:
                    continue
                terms, gene, _tx, biotype = p[0], p[1], p[2], p[3]
                seen_biotypes[biotype] += 1
                if biotype != "protein_coding":
                    dropped_biotype += 1
                    continue
                if gene not in panel:
                    n_nonpanel += 1
                    continue
                for t in terms.split("&"):
                    cls = classify(t)
                    if cls is None:
                        continue
                    cur = assign[vid].get(gene)
                    if cur is None or CLASS_RANK.index(cls) < CLASS_RANK.index(cur):
                        assign[vid][gene] = cls

    print("csq rows read              : %d" % n_rows)
    print("rows with a consequence    : %d" % n_with_csq)
    print("entries dropped, non-coding transcript biotype: %d" % dropped_biotype)
    print("entries dropped, gene not in panel            : %d" % n_nonpanel)
    print()
    print("transcript biotypes seen:")
    for b, c in sorted(seen_biotypes.items(), key=lambda x: -x[1])[:8]:
        print("  %-24s %d" % (b, c))
    print()

    # ------------------------------------------------------------- write files
    annot_path = "%s/regenie_annot.tsv" % args.outdir
    set_path = "%s/regenie_setlist.tsv" % args.outdir
    mask_path = "%s/regenie_masks.def" % args.outdir

    by_gene = defaultdict(list)
    n_annot = 0
    with open(annot_path, "w") as out:
        for vid in sorted(assign):
            for gene, cls in sorted(assign[vid].items()):
                out.write("%s %s %s\n" % (vid, gene, cls))
                by_gene[gene].append(vid)
                n_annot += 1

    with open(set_path, "w") as out:
        for gene in sorted(by_gene, key=lambda g: (len(anchor[g][0]), anchor[g])):
            chrom, pos = anchor[gene]
            vids = sorted(by_gene[gene],
                          key=lambda v: int(v.split(":")[1]))
            out.write("%s %s %d %s\n" % (gene, chrom, pos, ",".join(vids)))

    with open(mask_path, "w") as out:
        for name, cats in MASKS:
            out.write("%s %s\n" % (name, ",".join(cats)))

    print("wrote %s  (%d variant-gene lines)" % (annot_path, n_annot))
    print("wrote %s  (%d genes)" % (set_path, len(by_gene)))
    print("wrote %s  (%d masks)" % (mask_path, len(MASKS)))
    print()

    # -------------------------------------------------- testability by mask x bin
    counts = defaultdict(lambda: defaultdict(int))
    missing_maf = 0
    for vid, gm in assign.items():
        m = maf.get(vid)
        if m is None:
            missing_maf += 1
            continue
        for gene, cls in gm.items():
            for b in AAF_BINS:
                if m < b:
                    counts[(gene, b)][cls] += 1

    if missing_maf:
        print("variants with no frequency record: %d" % missing_maf)

    rows = []
    for (gene, b), c in counts.items():
        lof = c.get("LoF", 0)
        mis = c.get("missense", 0)
        syn = c.get("synonymous", 0)
        rows.append((gene, b, lof, mis, syn))

    with open("results/mask_counts.tsv", "w") as out:
        out.write("gene\taaf_bin\tM1_lof\tM2_lof_mis\tM4_all\n")
        for gene, b, lof, mis, syn in sorted(rows):
            out.write("%s\t%g\t%d\t%d\t%d\n" % (gene, b, lof, lof + mis,
                                                lof + mis + syn))

    print()
    print("=" * 66)
    print("TESTABLE GENES per mask x AAF bin")
    print("a gene needs at least 1 qualifying variant, and realistically more")
    print("=" * 66)
    print("%-10s %-12s %10s %10s %10s" % ("aaf bin", "min variants",
                                          "M1_lof", "M2_lof_mis", "M4_all"))
    print("-" * 66)
    for b in AAF_BINS:
        for floor in (1, 3, 10):
            n1 = sum(1 for g, bb, lof, mis, syn in rows if bb == b and lof >= floor)
            n2 = sum(1 for g, bb, lof, mis, syn in rows
                     if bb == b and lof + mis >= floor)
            n4 = sum(1 for g, bb, lof, mis, syn in rows
                     if bb == b and lof + mis + syn >= floor)
            print("%-10g %-12d %10d %10d %10d" % (b, floor, n1, n2, n4))
        print("-" * 66)

    print()
    print("genes with NO qualifying variant in any mask (check for annotation bugs):")
    got = {g for g, b, l, m, s in rows if l + m + s > 0}
    missing = sorted(panel - got)
    print("  %d: %s" % (len(missing), ", ".join(missing) if missing else "none"))
    print()
    print("wrote results/mask_counts.tsv")


if __name__ == "__main__":
    main()
