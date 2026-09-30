#!/usr/bin/env python3
"""
Day 1, step 1: turn the gene panel into GRCh38 regions, with locus strata.

Inputs
  config/panel_genes.txt          51 AD/PD shared modifier symbols, one per line
  config/panel_genes_17q21.txt    additional 17q21.31 candidates (not controls)
  ref/Homo_sapiens.GRCh38.<rel>.chr.gff3.gz

Outputs
  config/panel_genes_coords.tsv   symbol, coords, ENSG, biotype, list, cluster, kb
  config/panel_regions.bed        merged, padded, chr-prefixed slicing regions
  config/panel_regions_nochr.bed  same, Ensembl naming
  config/gene_strata.tsv          symbol -> stratum, for the analysis stage

Usage
  python3 scripts/01_make_regions.py --gff ref/Homo_sapiens.GRCh38.112.chr.gff3.gz
"""

import argparse
import gzip
import sys
from collections import defaultdict

MAIN_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

ALIASES = {
    "ERO1A": ["ERO1L"],
    "SLC49A3": ["MFSD7"],
    "HIGD1B": ["C17orf50"],
    "PLIN3": ["M6PRBP1", "TIP47"],
    "ARL17A": ["ARL17"],
    "LRRC37A2": ["LRRC37A2"],
    "STH": ["STH"],
}

# Gene clusters in this panel. Boundaries are set generously to cover the region
# of extended LD, they are not published coordinates. Genes inside a cluster are
# NOT independent tests, so results are stratified by cluster and confidence
# intervals are taken over replicates, never over genes. See docs/panel_notes.md.
CLUSTERS = [
    # name          chrom  start        end          slice_whole_interval
    ("17q21.31",    "17",  45_300_000,  46_900_000,  True),
    ("16p11.2",     "16",  29_600_000,  31_500_000,  False),
]


def parse_attrs(field):
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def load_genes(gff_path, wanted):
    lookup = defaultdict(list)
    search = set(wanted)
    for sym, alts in ALIASES.items():
        if sym in search:
            search.update(alts)

    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            if f[0] not in MAIN_CHROMS:
                continue
            attrs = parse_attrs(f[8])
            name = attrs.get("Name")
            if name is None or name not in search:
                continue
            lookup[name].append(
                (f[0], int(f[3]), int(f[4]),
                 attrs.get("ID", "").replace("gene:", ""),
                 attrs.get("biotype", "NA"))
            )
    return lookup


def pick_one(sym, hits):
    if len(hits) == 1:
        return hits[0], None
    coding = [h for h in hits if h[4] == "protein_coding"]
    pool = coding if coding else hits
    best = max(pool, key=lambda h: h[2] - h[1])
    note = "%s: %d GFF entries (%s), kept %s:%d-%d" % (
        sym, len(hits), ",".join("%s:%d" % (h[0], h[1]) for h in hits),
        best[0], best[1], best[2])
    return best, note


def cluster_of(chrom, start, end):
    for name, c, s, e, _ in CLUSTERS:
        if chrom == c and start < e and end > s:
            return name
    return "unclustered"


def merge(intervals, pad):
    by_chrom = defaultdict(list)
    for c, s, e in intervals:
        by_chrom[c].append((max(0, s - pad), e + pad))
    merged = []
    for c in sorted(by_chrom, key=lambda x: (len(x), x)):
        cur = None
        for s, e in sorted(by_chrom[c]):
            if cur and s <= cur[1]:
                cur = (cur[0], max(cur[1], e))
            else:
                if cur:
                    merged.append((c, cur[0], cur[1]))
                cur = (s, e)
        if cur:
            merged.append((c, cur[0], cur[1]))
    return merged


def read_list(path):
    out = []
    try:
        with open(path) as fh:
            for line in fh:
                line = line.split("#")[0].strip()
                if line:
                    out.append(line)
    except FileNotFoundError:
        pass
    return list(dict.fromkeys(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genes", default="config/panel_genes.txt")
    ap.add_argument("--extra", default="config/panel_genes_17q21.txt",
                    help="additional candidate genes; pass /dev/null to exclude")
    ap.add_argument("--pad", type=int, default=1000)
    ap.add_argument("--outdir", default="config")
    args = ap.parse_args()

    panel = read_list(args.genes)
    extra = read_list(args.extra)
    wanted = panel + [g for g in extra if g not in panel]
    membership = {g: "modifier" for g in panel}
    for g in extra:
        membership.setdefault(g, "candidate_17q21")

    print("shared modifier list    : %d" % len(panel))
    print("additional 17q21 candidates: %d" % len([g for g in extra if g not in panel]))
    print("panel total             : %d" % len(wanted))

    lookup = load_genes(args.gff, wanted)

    rows, missing, notes = [], [], []
    for sym in wanted:
        hits = lookup.get(sym)
        if not hits:
            for alt in ALIASES.get(sym, []):
                if lookup.get(alt):
                    hits = lookup[alt]
                    notes.append("%s matched via alias %s" % (sym, alt))
                    break
        if not hits:
            missing.append("%s (%s)" % (sym, membership[sym]))
            continue
        best, note = pick_one(sym, hits)
        if note:
            notes.append(note)
        rows.append((sym,) + best)

    coords_path = "%s/panel_genes_coords.tsv" % args.outdir
    strata_path = "%s/gene_strata.tsv" % args.outdir
    counts = defaultdict(int)
    with open(coords_path, "w") as out, open(strata_path, "w") as st:
        out.write("symbol\tchrom\tstart\tend\tensg\tbiotype\tlist\tcluster\tkb\n")
        st.write("symbol\tlist\tcluster\n")
        for sym, chrom, start, end, ensg, bt in sorted(rows, key=lambda r: (len(r[1]), r[1], r[2])):
            cl = cluster_of(chrom, start, end)
            ms = membership[sym]
            counts[(ms, cl)] += 1
            out.write("%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%.1f\n"
                      % (sym, chrom, start, end, ensg, bt, ms, cl, (end - start) / 1000))
            st.write("%s\t%s\t%s\n" % (sym, ms, cl))

    intervals = [(r[1], r[2], r[3]) for r in rows]
    for name, c, s, e, slice_whole in CLUSTERS:
        if slice_whole:
            intervals.append((c, s, e))
    merged = merge(intervals, args.pad)

    bed_chr = "%s/panel_regions.bed" % args.outdir
    bed_nochr = "%s/panel_regions_nochr.bed" % args.outdir
    with open(bed_chr, "w") as a, open(bed_nochr, "w") as b:
        for c, s, e in merged:
            a.write("chr%s\t%d\t%d\n" % (c, s, e))
            b.write("%s\t%d\t%d\n" % (c, s, e))

    total_bp = sum(e - s for _, s, e in merged)
    chroms = sorted({c for c, _, _ in merged}, key=lambda x: (len(x), x))
    kbs = sorted((r[3] - r[2]) / 1000 for r in rows)

    print()
    print("resolved : %d genes" % len(rows))
    print("MISSING  : %d %s" % (len(missing), missing if missing else ""))
    print()
    print("--- stratification (this is the number that matters) ---")
    for (ms, cl), n in sorted(counts.items()):
        print("  %-16s %-12s %3d genes" % (ms, cl, n))
    n_unclustered = sum(n for (_, cl), n in counts.items() if cl == "unclustered")
    n_clusters = len({cl for (_, cl) in counts if cl != "unclustered"})
    print("  effective independent units: ~%d (%d unclustered + %d clusters), not %d"
          % (n_unclustered + n_clusters, n_unclustered, n_clusters, len(rows)))
    print()
    print("--- gene size spread (power tracks variant count, which tracks size) ---")
    print("  min %.1f kb | median %.1f kb | max %.1f kb | ratio %.0fx"
          % (kbs[0], kbs[len(kbs) // 2], kbs[-1], kbs[-1] / kbs[0]))
    print()
    print("--- slicing ---")
    print("  %d merged intervals across %d chromosomes (%s)"
          % (len(merged), len(chroms), ",".join(chroms)))
    print("  total span %.2f Mb (pad %d bp each side)" % (total_bp / 1e6, args.pad))
    print()
    if notes:
        print("notes:")
        for n in notes:
            print("  " + n)
        print()
    for p in (coords_path, strata_path, bed_chr, bed_nochr):
        print("wrote %s" % p)

    if missing:
        print()
        print("ACTION: resolve missing symbols at https://www.genenames.org/")
        print("then add the Ensembl name to ALIASES in this script.")
        sys.exit(1)


if __name__ == "__main__":
    main()
