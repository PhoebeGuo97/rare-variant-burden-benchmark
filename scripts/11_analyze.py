#!/usr/bin/env python3
"""
Day 3, step 11: type I error and inflation from the null simulation grid.  v3.

Changes from v1:
  - refuses to read an arm that has both file naming schemes, instead of
    silently double counting the overlapping traits
  - Arm B is grouped by the REALIZED variance explained by superpopulation,
    not by the nominal strength label. The per-replicate population shift is
    drawn fresh, so nominal strength 0.50 produced realized R2 between 0.035
    and 0.386 across runs, and the three nominal groups overlap heavily.

Per mask x AAF bin x test it reports type I error at four alphas, genomic
inflation lambda, and the DISTRIBUTION of per-gene rates rather than only the
pooled rate. Genes here are not independent (14 in 16p11.2, 13 in 17q21.31,
about 31 independent units out of 56), so a pooled Wilson interval would assume
independence the panel does not have. Each gene's rate over independent
replicates is an honest binomial; those are computed first, then summarized.

Resolution: ~56 genes x N replicates per cell. At N=500, alpha 1e-3 is well
resolved and 1e-4 is marginal. Exome-wide 2.5e-6 is not estimable here and is
not reported as a point estimate.

No scipy: the chi-square quantile for df=1 is the squared normal quantile.

Usage
  python3 scripts/11_analyze_v2.py
"""

import argparse
import glob
import math
import gzip
import os
import re
from array import array
from collections import defaultdict
from statistics import NormalDist

ALPHAS = [0.05, 0.01, 1e-3, 1e-4]
CHI2_MEDIAN_1DF = 0.4549364231195736
TRAIT_RE = re.compile(r"^c\d+_")
ND = NormalDist()

ARMS = [
    ("a1_qt",           "A1  quantitative null, PCs"),
    ("a2_bt10_nofirth", "A2  binary 1:9,  no Firth"),
    ("a4_bt10_firth",   "A4  binary 1:9,  Firth"),
    ("a3_bt01_nofirth", "A3  binary 1:99, no Firth"),
    ("a5_bt01_firth",   "A5  binary 1:99, Firth"),
    ("b1_struct_nopc",  "B1  structured null, NO PCs"),
    ("b2_struct_pc",    "B2  structured null, PCs"),
]
STRUCT_ARMS = {"b1_struct_nopc", "b2_struct_pc"}


LN10 = math.log(10.0)
LN2 = math.log(2.0)
HALF_LN_2PI = 0.5 * math.log(2.0 * math.pi)


def chi2_from_log10p(lp):
    """Upper-tail chi-square with 1 df, from -log10(p). Never underflows.

    Going through p and calling inv_cdf(1 - p/2) breaks once p < ~2e-16,
    because 1 - p/2 rounds to exactly 1.0 in double precision. That is not a
    rare edge case here: the uncorrected 1:99 binary arm produces p-values far
    below that under a pure null, which is the whole point of running it.
    REGENIE reports LOG10P, so the conversion is done in log space instead.

    Below -log10(p) = 15 the exact normal quantile is used. Above it, solve
    q = phi(z)/z by fixed point. The two agree to ~0.1% at the crossover, and
    the asymptotic region only ever affects the extreme tail, never the median
    that lambda is computed from.
    """
    if lp != lp:
        return float("nan")
    if lp <= 0.0:
        return 0.0
    if lp < 15.0:
        return ND.inv_cdf(1.0 - (10.0 ** -lp) / 2.0) ** 2
    lq = -lp * LN10 - LN2
    z = math.sqrt(2.0 * -lq)
    for _ in range(80):
        zn = math.sqrt(max(2.0 * (-lq - math.log(z) - HALF_LN_2PI), 1e-12))
        if abs(zn - z) < 1e-13:
            z = zn
            break
        z = zn
    return z * z


def median_of(seq):
    s = sorted(x for x in seq if x == x)
    if not s:
        return float("nan")
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def quantiles(vals):
    s = sorted(v for v in vals if v == v)
    if not s:
        return (float("nan"),) * 3
    n = len(s)
    q = lambda f: s[min(n - 1, max(0, int(round(f * (n - 1)))))]
    return q(0.5), q(0.9), q(1.0)


def arm_files(tag, resdir):
    """Files for an arm, refusing to mix naming schemes.

    One run wrote <tag>_<TRAIT>.regenie.gz; the chunked rerun wrote
    <tag>_c<N>_<TRAIT>.regenie.gz. If a partial old run and a rerun both left
    files, a prefix glob reads the overlapping traits twice and overweights
    them, producing a plausible wrong number rather than a crash. Refuse.
    """
    files = sorted(glob.glob(os.path.join(resdir, "%s_*.regenie.gz" % tag)))
    chunked, plain = [], []
    for p in files:
        rest = os.path.basename(p)[len(tag) + 1:]
        (chunked if TRAIT_RE.match(rest) else plain).append(p)
    if chunked and plain:
        raise RuntimeError(
            "both naming schemes present: %d chunked, %d un-chunked. "
            "Overlapping traits would be counted twice. Delete the stale set."
            % (len(chunked), len(plain)))
    return chunked or plain


def trait_of(tag, path):
    rest = os.path.basename(path)[len(tag) + 1:].replace(".regenie.gz", "")
    return TRAIT_RE.sub("", rest)


def read_files(tag, files):
    for path in files:
        trait = trait_of(tag, path)
        with gzip.open(path, "rt") as fh:
            hdr = None
            for line in fh:
                if line.startswith("##"):
                    continue
                f = line.split()
                if hdr is None:
                    hdr = {k: i for i, k in enumerate(f)}
                    continue
                try:
                    lp = float(f[hdr["LOG10P"]])
                except (ValueError, KeyError, IndexError):
                    continue
                if lp != lp:
                    continue
                parts = f[hdr["ID"]].split(".")
                if len(parts) < 3:
                    continue
                # p may underflow to 0.0 for very large lp; that is fine for
                # the alpha comparisons, and lp carries the precision for chi2
                try:
                    p = 10.0 ** (-lp)
                except OverflowError:
                    p = 0.0
                yield (trait, parts[0], parts[1], ".".join(parts[2:]),
                       f[hdr["TEST"]], p, lp)


def realized_r2(pheno_path, ancestry_path):
    """Variance in each simulated trait explained by superpopulation."""
    anc = {}
    with open(ancestry_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                anc[f[0]] = f[2]
    with open(pheno_path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        traits = hdr[2:]
        sums = [defaultdict(float) for _ in traits]
        cnts = [defaultdict(int) for _ in traits]
        allv = [[] for _ in traits]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(hdr):
                continue
            pop = anc.get(f[1], "NA")
            for i, v in enumerate(f[2:]):
                try:
                    x = float(v)
                except ValueError:
                    continue
                sums[i][pop] += x
                cnts[i][pop] += 1
                allv[i].append(x)
    out = {}
    for i, t in enumerate(traits):
        vals = allv[i]
        if not vals:
            continue
        gm = sum(vals) / len(vals)
        sst = sum((v - gm) ** 2 for v in vals)
        ssb = sum(cnts[i][p] * (sums[i][p] / cnts[i][p] - gm) ** 2
                  for p in cnts[i] if cnts[i][p])
        out[t] = ssb / sst if sst else float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resdir", default="results/regenie")
    ap.add_argument("--struct-pheno", default="config/pheno_qt_struct_st.tsv")
    ap.add_argument("--ancestry", default="config/sample_ancestry.tsv")
    ap.add_argument("--out", default="results/calibration.tsv")
    ap.add_argument("--min-rep", type=int, default=20)
    args = ap.parse_args()

    group_of = lambda t: "-"
    if os.path.exists(args.struct_pheno) and os.path.exists(args.ancestry):
        r2 = realized_r2(args.struct_pheno, args.ancestry)
        vals = sorted(v for v in r2.values() if v == v)
        if vals:
            lo = vals[len(vals) // 3]
            hi = vals[2 * len(vals) // 3]
            print("Arm B realized R2 by superpopulation, across %d traits:" % len(vals))
            print("  min %.3f   tertile cuts %.3f / %.3f   max %.3f"
                  % (vals[0], lo, hi, vals[-1]))
            print("  grouping by realized R2, not by the nominal strength label,")
            print("  because the per-replicate shift is redrawn each time")
            print()

            def group_of(t, lo=lo, hi=hi, r2=r2):
                v = r2.get(t)
                if v is None or v != v:
                    return "R2_NA"
                return "R2_lo" if v < lo else ("R2_mid" if v < hi else "R2_hi")

    rows_out = []
    for tag, label in ARMS:
        try:
            files = arm_files(tag, args.resdir)
        except RuntimeError as e:
            print("SKIPPING %-20s %s" % (tag, e))
            continue
        if not files:
            print("%-22s no files" % tag)
            continue

        is_struct = tag in STRUCT_ARMS
        chis = defaultdict(lambda: array("d"))
        per_gene = defaultdict(lambda: [0] + [0] * len(ALPHAS))
        n_rows = 0
        for trait, gene, mask, abin, test, p, lp in read_files(tag, files):
            n_rows += 1
            cell = (mask, abin, test, group_of(trait) if is_struct else "-")
            chis[cell].append(chi2_from_log10p(lp))
            g = per_gene[cell + (gene,)]
            g[0] += 1
            for i, a in enumerate(ALPHAS):
                if p < a:
                    g[1 + i] += 1

        print("%-22s %4d files  %9d rows" % (tag, len(files), n_rows))

        genes_by_cell = defaultdict(list)
        for k, v in per_gene.items():
            if v[0] >= args.min_rep:
                genes_by_cell[k[:4]].append(v)

        for cell, arr in chis.items():
            mask, abin, test, grp = cell
            genes = genes_by_cell.get(cell, [])
            rec = {"arm": tag, "group": grp, "mask": mask, "aaf_bin": abin,
                   "test": test, "n_genes": len(genes), "n_tests": len(arr),
                   "lambda": (median_of(arr) / CHI2_MEDIAN_1DF
                              if len(arr) else float("nan"))}
            for i, a in enumerate(ALPHAS):
                num = sum(v[1 + i] for v in genes)
                den = sum(v[0] for v in genes)
                rec["t1e_%g" % a] = num / den if den else float("nan")
                med, p90, mx = quantiles([v[1 + i] / v[0] for v in genes if v[0]])
                rec["gene_med_%g" % a] = med
                rec["gene_p90_%g" % a] = p90
                rec["gene_max_%g" % a] = mx
            rows_out.append(rec)

    if not rows_out:
        print("nothing to analyze")
        return

    cols = (["arm", "group", "mask", "aaf_bin", "test", "n_genes", "n_tests",
             "lambda"]
            + ["t1e_%g" % a for a in ALPHAS]
            + ["gene_med_%g" % a for a in ALPHAS]
            + ["gene_p90_%g" % a for a in ALPHAS]
            + ["gene_max_%g" % a for a in ALPHAS])
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in rows_out:
            out.write("\t".join(
                ("%.5g" % r[c]) if isinstance(r.get(c), float) else str(r.get(c, ""))
                for c in cols) + "\n")
    print()
    print("wrote %s  (%d cells)" % (args.out, len(rows_out)))

    def show(tags, title, test="ADD", alpha=0.05, groups=None):
        rs = [r for r in rows_out if r["arm"] in tags and r["test"] == test
              and (groups is None or r["group"] in groups)]
        if not rs:
            return
        print()
        print("=" * 100)
        print("%s   [test=%s, alpha=%g]" % (title, test, alpha))
        print("=" * 100)
        print("%-20s %-7s %-12s %-10s %7s %8s %9s %9s %9s"
              % ("arm", "group", "mask", "aaf bin", "lambda", "T1E",
                 "gene med", "gene p90", "gene max"))
        print("-" * 100)
        for r in sorted(rs, key=lambda r: (r["arm"], r["group"], r["mask"],
                                           r["aaf_bin"])):
            print("%-20s %-7s %-12s %-10s %7.3f %8.4f %9.4f %9.4f %9.4f"
                  % (r["arm"], r["group"], r["mask"], r["aaf_bin"], r["lambda"],
                     r["t1e_%g" % alpha], r["gene_med_%g" % alpha],
                     r["gene_p90_%g" % alpha], r["gene_max_%g" % alpha]))

    print()
    print("Nominal type I error at alpha 0.05 is 0.05; lambda near 1.0 is calibrated.")
    print("'gene max' is the worst single gene. A mask can look fine pooled and be")
    print("badly anticonservative somewhere, which is the failure that matters.")

    show(["a1_qt"], "ARM A1: unstructured quantitative null, burden")
    show(["a2_bt10_nofirth", "a4_bt10_firth"], "BINARY 1:9, uncorrected vs Firth")
    show(["a3_bt01_nofirth", "a5_bt01_firth"], "BINARY 1:99, uncorrected vs Firth")
    show(["b1_struct_nopc", "b2_struct_pc"],
         "STRUCTURED NULL, without vs with PCs, most confounded tertile",
         groups={"R2_hi"})

    print()
    print("=" * 100)
    print("SKAT vs burden on the same null data (a1_qt).")
    print("Burden does not use the gamma code path that printed the boost warning,")
    print("so a SKAT-only miscalibration would point at that warning.")
    print("=" * 100)
    for test in ("ADD", "ADD-SKAT", "ADD-SKATO", "ADD-SKATO-ACAT",
                 "ADD-ACATV", "ADD-ACATO"):
        rs = [r for r in rows_out if r["arm"] == "a1_qt" and r["test"] == test]
        if not rs:
            continue
        print("  %-16s cells %3d  median lambda %6.3f  median T1E(0.05) %.4f"
              % (test, len(rs), median_of([r["lambda"] for r in rs]),
                 median_of([r["t1e_0.05"] for r in rs])))


if __name__ == "__main__":
    main()
