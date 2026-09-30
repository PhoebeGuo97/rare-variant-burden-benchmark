#!/usr/bin/env python3
"""
Day 2, step 6a: subset the Ensembl GFF3 to the panel genes and rename the
chromosomes to UCSC style.

Two things make this less trivial than a coordinate filter.

1. bcftools csq computes consequences from COMPLETE transcript models. A
   transcript whose CDS is clipped by a naive region filter would produce wrong
   consequences (a truncated CDS shifts the reading frame and invents stop
   codons). So this walks the GFF3 parent hierarchy: find genes overlapping the
   panel regions, then keep every descendant of those genes in full, wherever it
   lies.

2. Ensembl names chromosomes "17"; the 1000G VCF and our reference say "chr17".
   csq matches on contig name, so the GFF has to be renamed rather than the VCF,
   which would then disagree with the reference.

Input   ref/Homo_sapiens.GRCh38.<rel>.chr.gff3.gz, config/panel_regions_nochr.bed
Output  ref/panel.chr.gff3.gz

Usage
  python3 scripts/06_make_gff.py --gff ref/Homo_sapiens.GRCh38.112.chr.gff3.gz
"""

import argparse
import gzip
from collections import defaultdict


def parse_attrs(field):
    out = {}
    for kv in field.rstrip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def load_regions(path):
    reg = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                reg[f[0]].append((int(f[1]), int(f[2])))
    for c in reg:
        reg[c].sort()
    return reg


def overlaps(reg, chrom, start, end):
    for s, e in reg.get(chrom, []):
        if start <= e and end >= s:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--regions", default="config/panel_regions_nochr.bed")
    ap.add_argument("--out", default="ref/panel.chr.gff3")
    args = ap.parse_args()

    reg = load_regions(args.regions)
    print("regions: %d intervals on %d chromosomes"
          % (sum(len(v) for v in reg.values()), len(reg)))

    # pass 1: gene IDs overlapping the panel regions
    keep_genes = set()
    with gzip.open(args.gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            if overlaps(reg, f[0], int(f[3]), int(f[4])):
                gid = parse_attrs(f[8]).get("ID")
                if gid:
                    keep_genes.add(gid)
    print("genes overlapping panel regions: %d" % len(keep_genes))

    # pass 2: transcripts whose Parent is a kept gene
    keep_tx = set()
    with gzip.open(args.gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            parent = a.get("Parent", "")
            if parent in keep_genes and a.get("ID"):
                keep_tx.add(a["ID"])
    print("transcripts under those genes : %d" % len(keep_tx))

    # pass 3: write genes, their transcripts, and every child of those
    # transcripts, IN FULL, with the chr prefix added
    n_out = 0
    counts = defaultdict(int)
    with gzip.open(args.gff, "rt") as fh, open(args.out, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            fid, parent = a.get("ID", ""), a.get("Parent", "")
            if not (fid in keep_genes or fid in keep_tx
                    or parent in keep_genes or parent in keep_tx):
                continue
            f[0] = "chr" + f[0]
            out.write("\t".join(f) + "\n")
            counts[f[2]] += 1
            n_out += 1

    print()
    print("features written: %d" % n_out)
    for k in sorted(counts, key=lambda x: -counts[x])[:10]:
        print("  %-22s %d" % (k, counts[k]))
    print()
    if counts.get("CDS", 0) == 0:
        print("WARNING: no CDS features. csq will not produce protein")
        print("consequences. Check the GFF3 release and the region file.")
    print("wrote %s (bgzip next)" % args.out)


if __name__ == "__main__":
    main()
