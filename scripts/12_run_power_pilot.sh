#!/usr/bin/env bash
# Day 3, step 12: power grid PILOT.
#
# 80 traits, one gene, four effect sizes, both direction arms, 10 reps each.
# The only job is to find the h2 range where power is intermediate. A grid run
# at an h2 where everything is at 0% or everything is at 100% costs hours and
# shows nothing, because burden and SKAT cannot differ where both are saturated.
#
# Run:  conda activate rvb && bash scripts/12_run_power_pilot.sh

set -uo pipefail
ulimit -n 10240 2>/dev/null || true

P="data/plink/panel"
R="results/regenie"
COV="config/covariates.tsv"
PCS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
mkdir -p "$R" results data/plink

echo "=============== 1. export causal-gene genotypes ==============="
# Only the variants that can be causal: rare, coding, in the candidate genes.
awk '$3=="LoF" || $3=="missense" {print $1}' config/regenie_annot.tsv \
  | sort -u > data/plink/coding_ids.txt
echo "  rare+common coding variant IDs: $(wc -l < data/plink/coding_ids.txt | tr -d ' ')"

if [ ! -s data/plink/causal.raw ]; then
  plink2 --pfile "$P" \
         --extract data/plink/coding_ids.txt \
         --export A \
         --out data/plink/causal > data/plink/export.log 2>&1
  rc=$?
  if [ $rc -ne 0 ] || [ ! -s data/plink/causal.raw ]; then
    echo "  export FAILED:"; tail -15 data/plink/export.log; exit 1
  fi
fi
echo "  causal.raw: $(du -h data/plink/causal.raw | cut -f1)"

echo
echo "=============== 2. simulate ==============="
python3 scripts/12_sim_power.py --mode pilot || exit 1

echo
echo "=============== 3. REGENIE ==============="
t0=$(date +%s)
regenie \
  --step 2 --pgen "$P" \
  --phenoFile config/pheno_power_pilot.tsv \
  --covarFile "$COV" --covarColList "$PCS" \
  --anno-file config/regenie_annot.tsv \
  --set-list config/regenie_setlist.tsv \
  --mask-def config/regenie_masks.def \
  --aaf-bins 0.001,0.005,0.01 \
  --vc-tests skat,skato,skato-acat,acatv,acato-full \
  --build-mask max --minMAC 1 --ignore-pred --bsize 200 --gz \
  --out "$R/p0_pilot" > "$R/p0_pilot.stdout" 2>&1
echo "  exit $?  $(( $(date +%s) - t0 ))s  $(ls "$R"/p0_pilot_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ') files"
tail -4 "$R/p0_pilot.stdout"

echo
echo "=============== 4. where is power intermediate? ==============="
python3 scripts/13_power_analyze.py --tag p0_pilot --design results/power_design_pilot.tsv
