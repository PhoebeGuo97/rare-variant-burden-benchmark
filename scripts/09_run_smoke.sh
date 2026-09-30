#!/usr/bin/env bash
# Day 2, step 9: REGENIE smoke test.
#
# 20 phenotypes, not 1000. The point is to prove the file formats are accepted
# and to measure runtime before committing to the full grid. Generating 1000
# columns and discovering the set list has the wrong chromosome style would cost
# an hour.
#
# Run:  conda activate rvb && bash scripts/09_run_smoke.sh

set -uo pipefail

P="data/plink/panel"
mkdir -p results/regenie

echo "=============== 0. chromosome style check ==============="
# The set list names a chromosome per gene. REGENIE matches it against the pvar.
# plink2 may write "chr17" or "17" depending on its output-chr setting, and a
# mismatch silently yields zero tested sets rather than an error.
PVAR_CHR=$(grep -v '^#' "${P}.pvar" | head -1 | cut -f1)
SET_CHR=$(head -1 config/regenie_setlist.tsv | awk '{print $2}')
echo "  pvar chromosome : $PVAR_CHR"
echo "  set-list chrom  : $SET_CHR"

if [[ "$PVAR_CHR" == chr* && "$SET_CHR" != chr* ]]; then
  echo "  MISMATCH: adding chr prefix to the set list"
  awk '{$2="chr"$2; print}' config/regenie_setlist.tsv > config/regenie_setlist.fixed.tsv
  mv config/regenie_setlist.fixed.tsv config/regenie_setlist.tsv
elif [[ "$PVAR_CHR" != chr* && "$SET_CHR" == chr* ]]; then
  echo "  MISMATCH: stripping chr prefix from the set list"
  awk '{sub(/^chr/,"",$2); print}' config/regenie_setlist.tsv > config/regenie_setlist.fixed.tsv
  mv config/regenie_setlist.fixed.tsv config/regenie_setlist.tsv
else
  echo "  consistent"
fi

echo
echo "  variant ID match between pvar and annotation file:"
comm -12 <(grep -v '^#' "${P}.pvar" | cut -f3 | sort -u) \
         <(cut -d' ' -f1 config/regenie_annot.tsv | sort -u) | wc -l \
  | awk '{print "    IDs present in both: "$1}'
cut -d' ' -f1 config/regenie_annot.tsv | sort -u | wc -l \
  | awk '{print "    IDs in annotation file: "$1}'

echo
echo "=============== 1. phenotypes (20 reps) ==============="
python3 scripts/09_make_phenos.py --n-rep 20 --tag smoke || exit 1

echo
echo "=============== 2. REGENIE step 2, quantitative null ==============="
START=$(date +%s)
regenie \
  --step 2 \
  --pgen "$P" \
  --phenoFile config/pheno_qt_null_smoke.tsv \
  --covarFile config/covariates.tsv \
  --covarColList PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
  --anno-file config/regenie_annot.tsv \
  --set-list config/regenie_setlist.tsv \
  --mask-def config/regenie_masks.def \
  --aaf-bins 0.001,0.005,0.01 \
  --vc-tests skat,skato,skato-acat,acatv,acato-full \
  --build-mask max \
  --minMAC 1 \
  --ignore-pred \
  --bsize 200 \
  --gz \
  --out results/regenie/smoke_qt 2>&1 | tail -30
SEC=$(( $(date +%s) - START ))

echo
echo "  runtime: ${SEC}s for 20 quantitative traits"
echo "  extrapolated to 1000 traits: $(( SEC * 50 / 60 )) min (rough, scaling is not exactly linear)"

echo
echo "=============== 3. what came out ==============="
F=$(ls results/regenie/smoke_qt*.regenie.gz 2>/dev/null | head -1)
if [ -z "$F" ]; then
  echo "  NO OUTPUT. Check the tail above for the error."
  ls -la results/regenie/ 2>/dev/null
  exit 1
fi
echo "  file: $F"
echo "  header:"
zcat < "$F" | head -1 | tr ' ' '\n' | head -20 | awk '{print "    "$0}'
echo
echo "  rows: $(zcat < "$F" | tail -n +2 | wc -l | tr -d ' ')"
echo
echo "  tests present:"
zcat < "$F" | tail -n +2 | awk '{print $NF}' | sort | uniq -c | head -5
echo
echo "  mask x AAF combinations tested (from the ID column):"
zcat < "$F" | tail -n +2 | awk '{print $3}' | sed 's/^[^.]*\.//' | sort | uniq -c | sort -rn | head -20
echo
echo "DONE. If this looks right, the full run is the same command with"
echo "--n-rep 1000 and the three other phenotype files."
