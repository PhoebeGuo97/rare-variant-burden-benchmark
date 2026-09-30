#!/usr/bin/env bash
# Day 2, step 8: build the PLINK2 genotype set and the principal components.
#
# Run:  conda activate rvb && bash scripts/08_build_pgen.sh

set -uo pipefail

CSQ="data/csq"
OUT="data/plink"
mkdir -p "$OUT" results

CHROMS=$(cut -f1 config/panel_regions.bed | sort -u -V)

echo "=============== 1. concatenate ==============="
MERGED="$OUT/panel.vcf.gz"
if [ -s "${MERGED}.tbi" ]; then
  echo "  cached"
else
  LIST="$OUT/concat_list.txt"
  : > "$LIST"
  for c in $CHROMS; do echo "$CSQ/${c}.csq.vcf.gz" >> "$LIST"; done
  echo "  $(wc -l < "$LIST" | tr -d ' ') files, in chromosome order"
  bcftools concat -f "$LIST" -Oz -o "$MERGED" || exit 1
  tabix -f -p vcf "$MERGED"
fi
echo "  variants: $(bcftools index -n "$MERGED")"

echo
echo "=============== 2. pgen ==============="
# --double-id makes FID = IID = sample name, which is what the phenotype and
# covariate files will use. IDs are already chr:pos:ref:alt from step 6 and must
# NOT be regenerated here: REGENIE joins the annotation file, the set list and
# the genotypes by ID alone, and a second ID scheme silently drops every variant.
if [ -s "$OUT/panel.pgen" ]; then
  echo "  cached"
else
  plink2 --vcf "$MERGED" --double-id --vcf-half-call m \
         --make-pgen --sort-vars --out "$OUT/panel" || exit 1
fi
echo "  $(grep -vc '^#' "$OUT/panel.pvar") variants, $(grep -vc '^#' "$OUT/panel.psam") samples"

echo
echo "  ID scheme check (must be chr:pos:ref:alt):"
grep -v '^#' "$OUT/panel.pvar" | head -2 | awk '{print "    "$3}'
echo "  duplicate IDs: $(grep -v '^#' "$OUT/panel.pvar" | cut -f3 | sort | uniq -d | wc -l | tr -d ' ')"

echo
echo "=============== 3. PCA variant set ==============="
# PCs are computed from the panel slice itself, not from whole chromosomes.
# Downloading chr20-22 in full would be several GB for what is ultimately a few
# thousand pruned SNPs, and five superpopulations separate cleanly on far fewer.
#
# Two exclusions matter:
#   17q21.31, because a 1.5 Mb inversion polymorphism in strong LD would load
#   onto its own PC and describe an inversion rather than ancestry;
#   the seven VQSLOD < 0 genes, because their genotypes are the ones we just
#   showed are unreliable.
#
# Limitation to state in the write-up: these PCs come from 4.98 Mb of targeted
# gene regions, not a genome-wide scan. They are adequate for separating
# superpopulations, which is all the structure-confounded arm requires.
cat > "$OUT/pca_exclude.txt" <<'EOF'
PLEKHM1
ARL17B
ARL17A
LRRC37A
LRRC37A2
NSF
ZNF717
EOF

awk 'NR>1 && $8=="17q21.31" {print "chr"$2"\t"$3-1"\t"$4}' \
    config/panel_genes_coords.tsv > "$OUT/pca_exclude.bed"
awk 'NR==FNR{x[$1];next} FNR>1 && ($1 in x){print "chr"$2"\t"$3-1"\t"$4}' \
    "$OUT/pca_exclude.txt" config/panel_genes_coords.tsv >> "$OUT/pca_exclude.bed"
printf "chr17\t45300000\t46900000\n" >> "$OUT/pca_exclude.bed"
echo "  excluded intervals: $(wc -l < "$OUT/pca_exclude.bed" | tr -d ' ')"

plink2 --pfile "$OUT/panel" \
       --exclude range "$OUT/pca_exclude.bed" \
       --maf 0.05 --max-alleles 2 --snps-only \
       --indep-pairwise 200 50 0.2 \
       --out "$OUT/prune" >/dev/null 2>&1

echo "  common SNPs after pruning: $(wc -l < "$OUT/prune.prune.in" | tr -d ' ')"

echo
echo "=============== 4. PCA ==============="
plink2 --pfile "$OUT/panel" \
       --extract "$OUT/prune.prune.in" \
       --pca 10 --out "$OUT/pca" || exit 1

echo
echo "  variance explained (top 6):"
head -7 "$OUT/pca.eigenval" | awk '{printf "    PC%d  %.3f\n", NR, $1}'

echo
echo "=============== 5. covariate file ==============="
# FID IID PC1..PC10, tab separated, which is the format REGENIE expects.
{
  printf "FID\tIID"
  for i in $(seq 1 10); do printf "\tPC%d" "$i"; done
  printf "\n"
  tail -n +2 "$OUT/pca.eigenvec"
} > config/covariates.tsv
echo "  rows: $(( $(wc -l < config/covariates.tsv) - 1 ))"

echo
echo "  do the PCs separate the superpopulations? mean PC1 and PC2 by group:"
python3 - <<'PY'
import csv
from collections import defaultdict
anc={}
for line in open('config/sample_ancestry.tsv'):
    f=line.rstrip('\n').split('\t')
    if len(f)>=3: anc[f[0]]=f[2]
g=defaultdict(lambda: [[],[]])
with open('config/covariates.tsv') as fh:
    hdr=fh.readline().rstrip('\n').split('\t')
    i1,i2=hdr.index('PC1'),hdr.index('PC2')
    for line in fh:
        f=line.rstrip('\n').split('\t')
        s=anc.get(f[1])
        if s:
            g[s][0].append(float(f[i1])); g[s][1].append(float(f[i2]))
for k in sorted(g):
    a,b=g[k]
    print("    %-5s n=%4d  PC1 %+7.4f  PC2 %+7.4f" % (k,len(a),sum(a)/len(a),sum(b)/len(b)))
print()
print("    If the five groups have clearly distinct PC1/PC2 means, the PCs are")
print("    doing their job and the structure-confounded arm has something to correct.")
PY

echo
echo "DONE."
