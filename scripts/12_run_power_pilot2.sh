#!/usr/bin/env bash
# Day 3, step 12: power grid PILOT 2.
#
# Fixes two things pilot 1 exposed.
#
# 1. --vc-maxAAF 0.01
#    REGENIE computes the VC tests (SKAT, SKATO, ACAT) once per mask over
#    variants below --vc-maxAAF, and reports them in the "all" bin. The default
#    applies no AAF restriction, so pilot 1 was running SKAT on rare AND common
#    variants together, and no VC test appeared in any rare bin at all. Setting
#    it to 0.01 makes the SKAT-vs-burden comparison a rare variant comparison,
#    which is the only version of it worth reporting.
#
# 2. h2 grid moved down by roughly 10x
#    Every same-direction cell in pilot 1 sat at power 1.000, including the
#    lowest effect size. Two tests cannot be compared where both are saturated.
#
# Also tightens the mixed-direction sign balance; see BALANCE_TOL in the
# simulator.
#
# Run:  conda activate rvb && bash scripts/12_run_power_pilot2.sh

set -uo pipefail
ulimit -n 10240 2>/dev/null || true

P="data/plink/panel"
R="results/regenie"
COV="config/covariates.tsv"
PCS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
mkdir -p "$R" results

echo "=============== 1. simulate (lower h2, balanced mixed arm) ==============="
python3 scripts/12_sim_power.py --mode pilot2 || exit 1

echo
echo "=============== 2. REGENIE with --vc-maxAAF 0.01 ==============="
t0=$(date +%s)
regenie \
  --step 2 --pgen "$P" \
  --phenoFile config/pheno_power_pilot2.tsv \
  --covarFile "$COV" --covarColList "$PCS" \
  --anno-file config/regenie_annot.tsv \
  --set-list config/regenie_setlist.tsv \
  --mask-def config/regenie_masks.def \
  --aaf-bins 0.001,0.005,0.01 \
  --vc-tests skat,skato,skato-acat,acatv,acato-full \
  --vc-maxAAF 0.01 \
  --build-mask max --minMAC 1 --ignore-pred --bsize 200 --gz \
  --out "$R/p1_pilot2" > "$R/p1_pilot2.stdout" 2>&1
echo "  exit $?  $(( $(date +%s) - t0 ))s  $(ls "$R"/p1_pilot2_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ') files"
grep -i -E "vc-maxAAF|vc_maxAAF|maximum AAF|gene-based tests" "$R/p1_pilot2.stdout" | head -5
tail -3 "$R/p1_pilot2.stdout"

echo
echo "=============== 3. did VC tests appear on rare variants? ==============="
F=$(ls "$R"/p1_pilot2_*.regenie.gz 2>/dev/null | head -1)
if [ -n "$F" ]; then
  echo "  mask.bin -> tests present:"
  zcat < "$F" | tail -n +2 \
    | awk '{split($3,a,"."); key=a[2]; for(i=3;i<=length(a);i++) key=key"."a[i]; print key"\t"$8}' \
    | sort -u | awk -F'\t' '{m[$1]=m[$1]" "$2} END{for(k in m) printf "    %-22s%s\n", k, m[k]}' | sort
fi

echo
echo "=============== 4. power ==============="
python3 scripts/13_power_analyze.py --tag p1_pilot2 --design results/power_design_pilot2.tsv
