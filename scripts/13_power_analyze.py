#!/usr/bin/env python3
"""
Day 3, step 13: power, and the burden vs SKAT direction contrast.

Each trait has exactly one causal gene, so every result file splits cleanly:
  - the causal gene's rows measure POWER
  - the other 55 genes' rows measure type I error under a genetic alternative
    elsewhere in the panel, which is a free extra null check

Power is reported per mask x AAF bin x test, and indexed by the number of
qualifying variants, because that is the quantity the test consumes.

Usage
  python3 scripts/13_power_analyze.py --tag p0_pilot --design results/power_design_pilot.tsv
"""

import argparse
import glob
import gzip
import os
import re
from collections import defaultdict

TRAIT_RE = re.compile(r"^c\d+_")
ALPHAS = [0.05, 2.5e-6]


def load_design(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {k: n for n, k in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            d[f[i["trait"]]] = {
                "gene": f[i["causal_gene"]], "n_causal": int(f[i["n_causal"]]),
                "h2": float(f[i["h2"]]), "dir": f[i["direction"]],
            }
    return d


def trait_of(tag, path):
    rest = os.path.basename(path)[len(tag) + 1:].replace(".regenie.gz", "")
    return TRAIT_RE.sub("", rest)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--resdir", default="results/regenie")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    design = load_design(args.design)
    files = sorted(glob.glob(os.path.join(args.resdir, "%s_*.regenie.gz" % args.tag)))
    if not files:
        raise SystemExit("no files for tag %s" % args.tag)
    print("design traits %d, result files %d" % (len(design), len(files)))

    # (gene, h2, dir, mask, bin, test) -> [n, hits at each alpha]
    power = defaultdict(lambda: [0] + [0] * len(ALPHAS))
    # null rows from the non-causal genes
    null = defaultdict(lambda: [0] + [0] * len(ALPHAS))
    ncausal = {}

    for path in files:
        tr = trait_of(args.tag, path)
        d = design.get(tr)
        if d is None:
            continue
        ncausal[d["gene"]] = d["n_causal"]
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
                gene, mask, abin = parts[0], parts[1], ".".join(parts[2:])
                test = f[hdr["TEST"]]
                try:
                    p = 10.0 ** (-lp)
                except OverflowError:
                    p = 0.0
                if gene == d["gene"]:
                    k = (d["gene"], d["h2"], d["dir"], mask, abin, test)
                    tgt = power[k]
                else:
                    tgt = null[(d["dir"], mask, abin, test)]
                tgt[0] += 1
                for i, a in enumerate(ALPHAS):
                    if p < a:
                        tgt[1 + i] += 1

    out = args.out or "results/power_%s.tsv" % args.tag
    with open(out, "w") as fh:
        fh.write("gene\tn_causal\th2\tdirection\tmask\taaf_bin\ttest\tn\t"
                 + "\t".join("power_%g" % a for a in ALPHAS) + "\n")
        for (gene, h2, d, mask, abin, test), v in sorted(power.items()):
            fh.write("%s\t%d\t%g\t%s\t%s\t%s\t%s\t%d\t%s\n"
                     % (gene, ncausal.get(gene, 0), h2, d, mask, abin, test, v[0],
                        "\t".join("%.4f" % (v[1 + i] / v[0] if v[0] else float("nan"))
                                  for i in range(len(ALPHAS)))))
    print("wrote %s" % out)

    def rate(v, i=0):
        return v[1 + i] / v[0] if v[0] else float("nan")

    print()
    print("=" * 88)
    print("POWER at alpha 0.05, M2_lof_mis / 0.01, by effect size and direction")
    print("=" * 88)
    print("%-10s %6s %8s %12s %12s %10s" % ("gene", "nvar", "h2", "same-dir",
                                            "mixed-dir", "ratio"))
    print("-" * 88)
    seen = sorted({(g, h) for (g, h, d, m, b, t) in power
                   if m == "M2_lof_mis" and b == "0.01" and t == "ADD"})
    for g, h in seen:
        s = power.get((g, h, "same", "M2_lof_mis", "0.01", "ADD"))
        x = power.get((g, h, "mixed", "M2_lof_mis", "0.01", "ADD"))
        if not s or not x:
            continue
        ps, px = rate(s), rate(x)
        print("%-10s %6d %8g %12.3f %12.3f %10s"
              % (g, ncausal.get(g, 0), h, ps, px,
                 ("%.2f" % (ps / px)) if px > 0 else "inf"))

    print()
    print("=" * 88)
    print("THE CONTRAST: burden vs SKAT-O under mixed-direction effects")
    print("burden sums, so opposing signs cancel; SKAT tests variance and ignores sign")
    print("=" * 88)
    print("%-10s %8s %-7s %12s %12s %12s" % ("gene", "h2", "dir", "burden(ADD)",
                                             "SKAT", "SKATO"))
    print("-" * 88)
    for g, h in seen:
        for d in ("same", "mixed"):
            row = []
            for t in ("ADD", "ADD-SKAT", "ADD-SKATO"):
                v = power.get((g, h, d, "M2_lof_mis", "0.01", t))
                row.append(rate(v) if v else float("nan"))
            print("%-10s %8g %-7s %12.3f %12.3f %12.3f" % (g, h, d, *row))

    print()
    print("=" * 88)
    print("NULL CHECK: the 55 non-causal genes in the same runs, alpha 0.05")
    print("=" * 88)
    for d in ("same", "mixed"):
        v = null.get((d, "M2_lof_mis", "0.01", "ADD"))
        if v:
            print("  %-6s n=%7d  false positive rate %.4f  (nominal 0.05)"
                  % (d, v[0], rate(v)))

    print()
    print("If power is 0.00 or 1.00 everywhere, the h2 grid is in the wrong place")
    print("and burden vs SKAT cannot separate. Pick an h2 with power near 0.3-0.7")
    print("and run the full grid there.")


if __name__ == "__main__":
    main()
