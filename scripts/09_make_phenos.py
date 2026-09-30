#!/usr/bin/env python3
"""
Day 2, step 9: simulate phenotypes.

ARM A, unstructured null. y drawn independently of every genotype. Any rejection
is a false positive by construction. Quantitative, plus binary at two
case:control ratios, because unbalanced binary traits are where score tests
break and where Firth correction earns its place.

ARM B, structure confounded null. The phenotype mean is shifted by
superpopulation. Still null with respect to every gene, but now correlated with
ancestry, so any gene whose rare burden is concentrated in one superpopulation
will look associated unless structure is controlled.

  WHY NOT y = a*PC1 + noise, which was the original plan: correcting with the
  same PC1 used to build the phenotype makes the correction exact by
  construction, and the arm would demonstrate nothing. The superpopulation label
  is a different variable from the PCs, so whether the PCs absorb it is a real
  question with a real answer.

  What this arm is really testing: rare variants carry sharper population
  structure than common ones, and these PCs are built from common variants.
  So PC correction should be incomplete, and incomplete by MORE in the rarest
  AAF bin. If inflation after correction is worse at MAF < 0.001 than at
  MAF < 0.01, that is the result, and it is the reason structure is a bigger
  problem for rare variant tests than for GWAS.

Usage
  python3 scripts/09_make_phenos.py --n-rep 20 --tag smoke
  python3 scripts/09_make_phenos.py --n-rep 1000 --tag full
"""

import argparse
import random
from collections import defaultdict


def gauss_list(rng, n):
    return [rng.gauss(0.0, 1.0) for _ in range(n)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psam", default="data/plink/panel.psam")
    ap.add_argument("--ancestry", default="config/sample_ancestry.tsv")
    ap.add_argument("--n-rep", type=int, default=20)
    ap.add_argument("--tag", default="smoke")
    ap.add_argument("--seed", type=int, default=20260919)
    ap.add_argument("--outdir", default="config")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # sample order must follow the psam, and FID/IID must match it exactly
    samples = []
    with open(args.psam) as fh:
        for line in fh:
            if line.startswith("#"):
                head = line.lstrip("#").split()
                continue
            f = line.split()
            samples.append((f[0], f[1]))
    n = len(samples)
    print("samples: %d" % n)

    anc = {}
    with open(args.ancestry) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                anc[f[0]] = f[2]
    pops = sorted(set(anc.values()))
    print("superpopulations: %s" % ", ".join(pops))
    missing = sum(1 for _, iid in samples if iid not in anc)
    if missing:
        print("WARNING: %d samples have no ancestry label" % missing)

    R = args.n_rep
    tag = args.tag

    # ---------------------------------------------------- Arm A, quantitative
    path = "%s/pheno_qt_null_%s.tsv" % (args.outdir, tag)
    cols = ["QT%d" % (i + 1) for i in range(R)]
    data = [gauss_list(rng, n) for _ in range(R)]
    with open(path, "w") as out:
        out.write("FID\tIID\t" + "\t".join(cols) + "\n")
        for j, (fid, iid) in enumerate(samples):
            out.write("%s\t%s\t%s\n"
                      % (fid, iid, "\t".join("%.6f" % data[i][j] for i in range(R))))
    print("wrote %s  (%d traits)" % (path, R))

    # ---------------------------------------------------- Arm A, binary
    # 1:9 and 1:99. The second is the one that breaks uncorrected score tests.
    for prev, name in ((0.10, "bt10"), (0.01, "bt01")):
        path = "%s/pheno_%s_null_%s.tsv" % (args.outdir, name, tag)
        cols = ["%s_%d" % (name.upper(), i + 1) for i in range(R)]
        # fixed case count per replicate, so the ratio is exactly as advertised
        ncase = max(1, int(round(prev * n)))
        with open(path, "w") as out:
            out.write("FID\tIID\t" + "\t".join(cols) + "\n")
            mat = []
            for _ in range(R):
                y = [0] * n
                for idx in rng.sample(range(n), ncase):
                    y[idx] = 1
                mat.append(y)
            for j, (fid, iid) in enumerate(samples):
                out.write("%s\t%s\t%s\n"
                          % (fid, iid, "\t".join(str(mat[i][j]) for i in range(R))))
        print("wrote %s  (%d traits, %d cases / %d controls)"
              % (path, R, ncase, n - ncase))

    # ---------------------------------------------------- Arm B, structured
    # y = superpopulation shift + noise, at three confounding strengths.
    # The shift per population is drawn once per replicate so the direction is
    # not always the same population being high.
    path = "%s/pheno_qt_struct_%s.tsv" % (args.outdir, tag)
    strengths = [0.25, 0.50, 1.00]
    cols, mats = [], []
    for s in strengths:
        for i in range(R):
            shift = {p: rng.gauss(0.0, 1.0) for p in pops}
            col = []
            for _, iid in samples:
                col.append(s * shift.get(anc.get(iid, pops[0]), 0.0) + rng.gauss(0, 1))
            cols.append("STR%03d_%d" % (int(s * 100), i + 1))
            mats.append(col)
    with open(path, "w") as out:
        out.write("FID\tIID\t" + "\t".join(cols) + "\n")
        for j, (fid, iid) in enumerate(samples):
            out.write("%s\t%s\t%s\n"
                      % (fid, iid, "\t".join("%.6f" % m[j] for m in mats)))
    print("wrote %s  (%d traits: %d strengths x %d reps)"
          % (path, len(cols), len(strengths), R))

    # sanity: how much variance does the population label explain at each strength
    print()
    print("Arm B check, variance explained by superpopulation (first rep each):")
    for k, s in enumerate(strengths):
        col = mats[k * R]
        grp = defaultdict(list)
        for (fid, iid), v in zip(samples, col):
            grp[anc.get(iid, "NA")].append(v)
        gm = sum(col) / n
        ssb = sum(len(v) * (sum(v) / len(v) - gm) ** 2 for v in grp.values())
        sst = sum((v - gm) ** 2 for v in col)
        print("  strength %.2f  R2 = %.3f" % (s, ssb / sst))
    print()
    print("Without PC covariates these should inflate. With them, the question is")
    print("whether the inflation clears, and whether it clears equally in every")
    print("AAF bin. Expect the rarest bin to clear least.")


if __name__ == "__main__":
    main()
