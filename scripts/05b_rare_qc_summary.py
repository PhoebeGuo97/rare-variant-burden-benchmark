#!/usr/bin/env python3
"""
Day 2, step 5b: per-gene QC computed on the RARE variants that will actually
enter the masks, and the corrected test of the mappability prediction.

The earlier version tested the prediction with heterozygote excess at MAF >= 0.05.
That set is disjoint from the variants under test (masks go to MAF < 0.01), and
het excess has no power at rare frequencies anyway, because it is detected
through missing homozygotes and the expected homozygote count is N*p^2.

This version restricts to MAF < --max-maf and uses frequency independent
metrics: mapping quality, mapping-quality-zero reads, Mendelian error rate from
the 1000G trios, and VQSLOD.

It also reports the quantity a burden test literally sums: rare alternate
alleles per sample per kb. If mismapping inflates that in a gene, the burden
test for that gene is broken in the way that matters, regardless of what any
other metric says.

Usage
  python3 scripts/05b_rare_qc_summary.py
  python3 scripts/05b_rare_qc_summary.py --max-maf 0.001
"""

import argparse
import gzip
import sys
from collections import defaultdict

N_SAMPLES = 2504


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                yield dict(zip(header, f))


def num(x):
    if x is None or x == "." or x == "":
        return None
    try:
        v = float(x)
        return None if v != v else v
    except ValueError:
        return None


def median(xs):
    xs = sorted(x for x in xs if x is not None and x == x)
    if not xs:
        return float("nan")
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def pct(n, d):
    return 100.0 * n / d if d else float("nan")


def load_genes(path):
    genes, meta = defaultdict(list), {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = {k: n for n, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            sym = f[i["symbol"]]
            genes["chr" + f[i["chrom"]]].append(
                (int(f[i["start"]]), int(f[i["end"]]), sym))
            meta[sym] = {"cluster": f[i["cluster"]], "list": f[i["list"]],
                         "kb": float(f[i["kb"]])}
    for c in genes:
        genes[c].sort()
    return genes, meta


def load_map(path):
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i = {k: n for n, k in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= len(header):
                out[f[i["symbol"]]] = (float(f[i["segdup_frac"]]),
                                       f[i["segdup_flag"]])
    return out


def genes_at(genes, chrom, pos):
    return [s for a, b, s in genes.get(chrom, []) if a <= pos <= b]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", default="config/panel_genes_coords.tsv")
    ap.add_argument("--map", default="config/gene_mappability.tsv")
    ap.add_argument("--freq", default="results/variant_qc_pooled.tsv.gz",
                    help="from scripts/04_variant_qc.sh, supplies AC/AN/F_MISSING")
    ap.add_argument("--annot", default="results/site_annot.tsv.gz",
                    help="from scripts/05_rare_qc.sh")
    ap.add_argument("--max-maf", type=float, default=0.01,
                    help="upper MAF for the rare set under test (default 0.01)")
    ap.add_argument("--out", default="results/gene_rare_qc.tsv")
    args = ap.parse_args()

    genes, meta = load_genes(args.coords)
    mapp = load_map(args.map)

    try:
        freq = {}
        for r in read_tsv_gz(args.freq):
            ac, an = num(r["AC"]), num(r["AN"])
            freq[(r["CHROM"], r["POS"], r["REF"], r["ALT"])] = (
                ac, an, num(r["F_MISSING"]))
    except FileNotFoundError:
        sys.exit("missing %s, run scripts/04_variant_qc.sh first" % args.freq)

    agg = defaultdict(lambda: {
        "n_rare": 0, "ac_sum": 0.0, "mq": [], "mq0_pos": 0, "me": [],
        "me_pos": 0, "me_n": 0, "vqslod": [], "fmiss": [], "mq_low": 0,
        "inb": []})

    n_annot = n_rare_total = n_nofreq = 0
    for r in read_tsv_gz(args.annot):
        n_annot += 1
        key = (r["CHROM"], r["POS"], r["REF"], r["ALT"])
        f = freq.get(key)
        if f is None:
            n_nofreq += 1
            continue
        ac, an, fmiss = f
        if ac is None or an is None or an <= 0:
            continue
        af = ac / an
        maf = min(af, 1.0 - af)
        if maf >= args.max_maf or ac < 1:
            continue
        n_rare_total += 1

        syms = genes_at(genes, r["CHROM"], int(r["POS"]))
        if not syms:
            continue
        mq, mq0 = num(r["MQ"]), num(r["MQ0"])
        me, vq = num(r["ME"]), num(r["VQSLOD"])
        inb = num(r["InbreedingCoeff"])
        for sym in syms:
            a = agg[sym]
            a["n_rare"] += 1
            a["ac_sum"] += ac
            if fmiss is not None:
                a["fmiss"].append(fmiss)
            if mq is not None:
                a["mq"].append(mq)
                if mq < 50:
                    a["mq_low"] += 1
            if mq0 is not None and mq0 > 0:
                a["mq0_pos"] += 1
            if vq is not None:
                a["vqslod"].append(vq)
            if inb is not None:
                a["inb"].append(inb)
            if me is not None:
                a["me_n"] += 1
                a["me"].append(me)
                if me > 0:
                    a["me_pos"] += 1

    print("site annotations read : %d" % n_annot)
    print("no frequency record    : %d" % n_nofreq)
    print("rare set (MAF < %.3f)  : %d variants" % (args.max_maf, n_rare_total))
    print()

    rows = []
    for sym, a in agg.items():
        frac, flag = mapp.get(sym, (float("nan"), "NA"))
        kb = meta[sym]["kb"]
        rows.append({
            "symbol": sym, "list": meta[sym]["list"],
            "cluster": meta[sym]["cluster"], "segdup_frac": frac,
            "segdup_flag": flag, "kb": kb,
            "n_rare": a["n_rare"],
            "rare_per_kb": a["n_rare"] / kb if kb else float("nan"),
            # the quantity a burden test sums
            "rare_alleles_per_sample_per_kb":
                a["ac_sum"] / N_SAMPLES / kb if kb else float("nan"),
            "median_MQ": median(a["mq"]),
            "pct_MQ_lt50": pct(a["mq_low"], a["n_rare"]),
            "pct_MQ0_gt0": pct(a["mq0_pos"], a["n_rare"]),
            "n_ME": a["me_n"],
            "pct_ME_gt0": pct(a["me_pos"], a["me_n"]),
            "median_ME": median(a["me"]),
            "median_VQSLOD": median(a["vqslod"]),
            "median_InbCoeff": median(a["inb"]),
            "median_fmiss": median(a["fmiss"]),
        })

    cols = ["symbol", "list", "cluster", "segdup_frac", "segdup_flag", "kb",
            "n_rare", "rare_per_kb", "rare_alleles_per_sample_per_kb",
            "median_MQ", "pct_MQ_lt50", "pct_MQ0_gt0", "n_ME", "pct_ME_gt0",
            "median_ME", "median_VQSLOD", "median_InbCoeff", "median_fmiss"]
    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in sorted(rows, key=lambda x: -x["segdup_frac"]):
            out.write("\t".join(("%.4f" % r[c]) if isinstance(r[c], float)
                                else str(r[c]) for c in cols) + "\n")

    print("=" * 96)
    print("PREDICTION TEST on the rare set (MAF < %.3f), the variants actually"
          " under test" % args.max_maf)
    print("=" * 96)
    hdr = ("%-13s %-11s %6s %7s %8s %8s %8s %8s %9s %9s"
           % ("gene", "cluster", "segdup", "n_rare", "medMQ", "MQ<50%",
              "MQ0>0%", "ME>0%", "medVQSLOD", "rare/kb"))
    print(hdr)
    print("-" * 96)

    def line(r):
        return ("%-13s %-11s %6.3f %7d %8.1f %8.1f %8.1f %8.1f %9.2f %9.2f"
                % (r["symbol"], r["cluster"], r["segdup_frac"], r["n_rare"],
                   r["median_MQ"], r["pct_MQ_lt50"], r["pct_MQ0_gt0"],
                   r["pct_ME_gt0"], r["median_VQSLOD"], r["rare_per_kb"]))

    for r in sorted(rows, key=lambda x: (-x["segdup_frac"], x["symbol"]))[:12]:
        print(line(r))
    print("-" * 96)
    print("predicted to look ORDINARY:")
    for sym in ("MAPT", "CRHR1", "KANSL1", "WNT3", "FUS", "ITGAM", "INPP5D"):
        r = next((x for x in rows if x["symbol"] == sym), None)
        if r:
            print(line(r))
    print()

    hi = [r for r in rows if r["segdup_flag"] == "high" and r["n_rare"] >= 20]
    lo = [r for r in rows if r["segdup_flag"] == "low" and r["n_rare"] >= 20]
    for label, grp in (("segdup high", hi), ("segdup low", lo)):
        if grp:
            print("%-12s n=%2d | medMQ %5.1f | MQ<50%% %5.1f | MQ0>0%% %5.1f | "
                  "ME>0%% %5.1f | rare alleles/sample/kb %6.3f"
                  % (label, len(grp),
                     median([g["median_MQ"] for g in grp]),
                     median([g["pct_MQ_lt50"] for g in grp]),
                     median([g["pct_MQ0_gt0"] for g in grp]),
                     median([g["pct_ME_gt0"] for g in grp]),
                     median([g["rare_alleles_per_sample_per_kb"] for g in grp])))
    print()
    if hi and lo:
        checks = []
        mq_h, mq_l = (median([g["median_MQ"] for g in hi]),
                      median([g["median_MQ"] for g in lo]))
        checks.append(("mapping quality lower in flagged genes", mq_h < mq_l - 2))
        m0_h, m0_l = (median([g["pct_MQ0_gt0"] for g in hi]),
                      median([g["pct_MQ0_gt0"] for g in lo]))
        checks.append(("more MQ0 reads in flagged genes", m0_h > m0_l + 5))
        me_h, me_l = (median([g["pct_ME_gt0"] for g in hi]),
                      median([g["pct_ME_gt0"] for g in lo]))
        if me_h == me_h and me_l == me_l:
            checks.append(("more Mendelian errors in flagged genes", me_h > me_l * 2))
        bu_h, bu_l = (median([g["rare_alleles_per_sample_per_kb"] for g in hi]),
                      median([g["rare_alleles_per_sample_per_kb"] for g in lo]))
        checks.append(("inflated rare allele burden in flagged genes", bu_h > bu_l * 1.5))

        print("VERDICT, criterion by criterion:")
        for name, ok in checks:
            print("  [%s] %s" % ("PASS" if ok else "fail", name))
        n_pass = sum(1 for _, ok in checks if ok)
        print()
        if n_pass >= 3:
            print("Prediction HELD on %d of %d criteria, measured on the rare"
                  " variants themselves." % (n_pass, len(checks)))
            print("The segdup flag is a usable stratum and no longer rests on a")
            print("common-variant proxy.")
        elif n_pass >= 1:
            print("MIXED: %d of %d. Report which metrics separated and which did"
                  " not." % (n_pass, len(checks)))
            print("Prefer the measured columns over the segdup flag for the stratum.")
        else:
            print("Prediction FAILED on the rare set. The segdup flag does not")
            print("predict rare variant quality here. Drop it as a stratum and say so.")
    print()
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
