#!/usr/bin/env python3
"""
Day 3, step 12: simulate phenotypes with real genetic effects, for the power grid.

The contrast this exists to produce: burden collapses a gene to a single sum,
so effects of opposite sign cancel and burden loses its power entirely. SKAT
tests a variance component and does not care about sign. That is the actual
answer to "how does SKAT-O differ from burden", and it is worth measuring rather
than asserting.

MODEL
  y_i = c * sum_j s_j g_ij + e_i,  e ~ N(0,1)

  s_j is +1 for every causal variant (SAME direction) or +/-1 balanced to
  near zero net sign (MIXED direction; see BALANCE_TOL). c is chosen per replicate so that the genetic component
  has variance h2 in every arm.

  Fixing variance explained rather than fixing the per-allele effect is what
  makes the comparison fair. With mixed signs the burden score cancels, so if
  you fixed beta instead, the mixed arm would carry less signal and burden
  would lose for two reasons at once, confounding "burden cannot see cancelling
  effects" with "there was less to see".

DESIGN
  Each phenotype column has exactly ONE causal gene. The other 55 genes in the
  same REGENIE run are null by construction, so every replicate yields one power
  observation and 55 simultaneous null observations at no extra cost.

Causal genes are chosen unclustered and with median VQSLOD > 0, spanning
qualifying variant counts, so that power can be indexed by variant count rather
than by gene.

Usage
  python3 scripts/12_sim_power.py --mode pilot
  python3 scripts/12_sim_power.py --mode full
"""

import argparse
import gzip
import os
import random
from collections import defaultdict

# unclustered, VQSLOD > 0, spanning rare coding variant counts at MAF < 0.01
CAUSAL_GENES = ["ABCE1", "ERO1A", "DLGAP1", "FBN1"]

# Pilot 1 put every same-direction cell at power 1.000, so the grid moved down.
MODES = {
    "pilot":  {"genes": ["DLGAP1"], "h2": [0.005, 0.01, 0.02, 0.04],
               "dirs": ["same", "mixed"], "reps": 10},
    "pilot2": {"genes": ["DLGAP1"], "h2": [0.0005, 0.001, 0.002, 0.004],
               "dirs": ["same", "mixed"], "reps": 15},
    # pilot3 sits where pilot2 showed power is intermediate, with more reps,
    # so the MAC-collapse comparison is not read off 15-replicate noise
    "pilot3": {"genes": ["DLGAP1"], "h2": [0.002, 0.004],
               "dirs": ["same", "mixed"], "reps": 25},
    "full":   {"genes": CAUSAL_GENES, "h2": [0.00125, 0.0025, 0.005],
               "dirs": ["same", "mixed"], "reps": 50},
}

# How close to balanced a "mixed" sign vector has to be, as a fraction of the
# number of causal variants. Pilot 1 used 0.5, which with 40 variants allows a
# 30/10 split; the residual net direction left burden with 0.90 power at high
# h2, so the arm was not testing what it claimed. At 0.05 the sum is near zero
# and burden has nothing to collapse onto.
#
# Note this makes "mixed" an idealized case. Truly random signs would sit
# between the two arms, because a random draw usually has some net imbalance.
BALANCE_TOL = 0.05


def load_annot(path, keep_cats=("LoF", "missense")):
    gene_vars = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3 and f[2] in keep_cats:
                gene_vars[f[1]].append(f[0])
    return gene_vars


def load_maf(path):
    maf = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            try:
                ac, an = float(f[4]), float(f[5])
            except ValueError:
                continue
            if an <= 0:
                continue
            af = ac / an
            maf["%s:%s:%s:%s" % (f[0], f[1], f[2], f[3])] = min(af, 1 - af)
    return maf


def read_raw(path):
    """plink2 --export A output: header then FID IID PAT MAT SEX PHENO g1 g2 ...

    Missing genotypes are written as NA; they are mean-imputed here, which is
    what REGENIE does internally for mask building, so the simulated effect and
    the tested genotype see the same values.
    """
    with open(path) as fh:
        hdr = fh.readline().split()
        varcols = hdr[6:]
        # plink appends _<allele> to each variant ID in the .raw header
        names = [c.rsplit("_", 1)[0] for c in varcols]
        ids, rows = [], []
        for line in fh:
            f = line.split()
            if len(f) < 6:
                continue
            ids.append((f[0], f[1]))
            rows.append([None if x == "NA" else float(x) for x in f[6:]])
    n = len(rows)
    m = len(names)
    for j in range(m):
        vals = [rows[i][j] for i in range(n) if rows[i][j] is not None]
        mu = sum(vals) / len(vals) if vals else 0.0
        for i in range(n):
            if rows[i][j] is None:
                rows[i][j] = mu
    return ids, names, rows


def variance(xs):
    n = len(xs)
    if n < 2:
        return 0.0
    m = sum(xs) / n
    return sum((x - m) ** 2 for x in xs) / (n - 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=sorted(MODES), default="pilot")
    ap.add_argument("--raw", default="data/plink/causal.raw")
    ap.add_argument("--annot", default="config/regenie_annot.tsv")
    ap.add_argument("--maf", default="results/variant_qc_pooled.tsv.gz")
    ap.add_argument("--max-maf", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=20260920)
    ap.add_argument("--outdir", default="config")
    args = ap.parse_args()

    cfg = MODES[args.mode]
    rng = random.Random(args.seed)

    gene_vars = load_annot(args.annot)
    maf = load_maf(args.maf)
    ids, names, G = read_raw(args.raw)
    col = {nm: j for j, nm in enumerate(names)}
    n = len(ids)
    print("samples %d, genotype columns %d" % (n, len(names)))

    design = []
    cols = []
    for gene in cfg["genes"]:
        vs = [v for v in gene_vars.get(gene, [])
              if v in col and maf.get(v, 1.0) < args.max_maf]
        if not vs:
            print("  %-10s no qualifying variants, skipped" % gene)
            continue
        idx = [col[v] for v in vs]
        print("  %-10s %d causal variants at MAF < %g" % (gene, len(vs), args.max_maf))
        for h2 in cfg["h2"]:
            for d in cfg["dirs"]:
                for r in range(cfg["reps"]):
                    s = ([1.0] * len(idx) if d == "same"
                         else [rng.choice((1.0, -1.0)) for _ in idx])
                    # mixed with an even split is the interesting case; reroll
                    # a degenerate all-same draw so the arm means what it says
                    if d == "mixed" and len(idx) > 3:
                        tol = max(2.0, BALANCE_TOL * len(idx))
                        tries = 0
                        while abs(sum(s)) > tol and tries < 10000:
                            s = [rng.choice((1.0, -1.0)) for _ in idx]
                            tries += 1
                    gv = [sum(G[i][j] * sj for j, sj in zip(idx, s)) for i in range(n)]
                    vg = variance(gv)
                    if vg <= 0:
                        continue
                    c = (h2 / vg) ** 0.5
                    y = [c * gv[i] + rng.gauss(0.0, 1.0) for i in range(n)]
                    name = "P_%s_h%g_%s_%d" % (gene, h2, d, r + 1)
                    cols.append((name, y))
                    design.append((name, gene, len(idx), h2, d, r + 1,
                                   sum(1 for x in s if x > 0),
                                   sum(1 for x in s if x < 0)))

    if not cols:
        raise SystemExit("no phenotypes generated")

    path = "%s/pheno_power_%s.tsv" % (args.outdir, args.mode)
    with open(path, "w") as out:
        out.write("FID\tIID\t" + "\t".join(c[0] for c in cols) + "\n")
        for i, (fid, iid) in enumerate(ids):
            out.write("%s\t%s\t%s\n"
                      % (fid, iid, "\t".join("%.6f" % c[1][i] for c in cols)))
    print("wrote %s  (%d traits)" % (path, len(cols)))

    dpath = "results/power_design_%s.tsv" % args.mode
    os.makedirs("results", exist_ok=True)
    with open(dpath, "w") as out:
        out.write("trait\tcausal_gene\tn_causal\th2\tdirection\trep\tn_pos\tn_neg\n")
        for row in design:
            out.write("\t".join(str(x) for x in row) + "\n")
    print("wrote %s" % dpath)

    print()
    print("Every trait has ONE causal gene; the other 55 genes are null in the")
    print("same run, so each replicate also contributes null observations.")


if __name__ == "__main__":
    main()
