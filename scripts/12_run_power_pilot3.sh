#!/usr/bin/env bash
# Day 3, step 12: power grid PILOT 3.
#
# Pilot 2 showed SKAT returning IDENTICAL statistics to burden for M1_lof
# (LOG10P 0.543054 and CHISQ 1.13655 for all six tests). The REGENIE log gives
# the reason:
#
#     -variants with MAC <= 10 are collapsed into a mask
#
# That is --vc-MACthr, default 10. Every variant with minor allele count at or
# below 10 is collapsed into one burden pseudo-variant BEFORE any
# variance-component test runs. At n = 2504, MAC 10 is MAF 0.002, and about 86%
# of this panel's rare variants fall below it. SKAT is therefore being handed a
# collapsed burden term plus a few individual variants, and the collapse is
# exactly what destroys the sign information SKAT exists to use.
#
# This script runs the SAME phenotypes under both settings so the default is
# measured rather than described:
#
#   arm A: --vc-MACthr 10   (REGENIE default)
#   arm B: --vc-MACthr 1    (no collapsing)
#
# The expected contrast, if this diagnosis is right: under MACthr 1, SKAT
# recovers power in the mixed-direction arm while burden does not. Under
# MACthr 10, the two stay locked together.
#
# If SKAT does NOT separate even at MACthr 1, the diagnosis is wrong and the
# collapse was not the explanation. That outcome gets reported too.
#
# Run:  conda activate rvb && bash scripts/12_run_power_pilot3.sh

set -uo pipefail
ulimit -n 10240 2>/dev/null || true

P="data/plink/panel"
R="results/regenie"
COV="config/covariates.tsv"
PCS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
mkdir -p "$R" results

echo "=============== 1. simulate ==============="
python3 scripts/12_sim_power.py --mode pilot3 || exit 1

run_one() {   # run_one <tag> <macthr>
  local tag=$1 mac=$2
  echo
  echo "=============== REGENIE, --vc-MACthr $mac ==============="
  local t0=$(date +%s)
  regenie \
    --step 2 --pgen "$P" \
    --phenoFile config/pheno_power_pilot3.tsv \
    --covarFile "$COV" --covarColList "$PCS" \
    --anno-file config/regenie_annot.tsv \
    --set-list config/regenie_setlist.tsv \
    --mask-def config/regenie_masks.def \
    --aaf-bins 0.001,0.005,0.01 \
    --vc-tests skat,skato,skato-acat,acatv,acato-full \
    --vc-maxAAF 0.01 \
    --vc-MACthr "$mac" \
    --build-mask max --minMAC 1 --ignore-pred --bsize 200 --gz \
    --out "$R/$tag" > "$R/${tag}.stdout" 2>&1
  echo "  exit $?  $(( $(date +%s) - t0 ))s  $(ls "$R/${tag}"_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ') files"
  grep -i "collapsed into a mask" "$R/${tag}.stdout" | head -1 | sed 's/^/  /'
}

run_one p2_mac10 10
run_one p2_mac1 1

echo
echo "=============== are SKAT and burden still identical? ==============="
# The decisive check: count DLGAP1 rows where SKAT and burden give the same
# LOG10P to 6 significant figures. If the collapse is the cause, this should be
# common at MACthr 10 and rare at MACthr 1.
for tag in p2_mac10 p2_mac1; do
  F=$(ls "$R/${tag}"_*.regenie.gz 2>/dev/null | head -1)
  [ -z "$F" ] && continue
  echo "  $tag:"
  zcat < "$F" | awk 'NR>1 && $3 ~ /^DLGAP1\./ {key=$3"|"$8; p[$3"|"$8]=$12}
    END{ for(k in p){ split(k,a,"|"); v[a[1]][a[2]]=p[k] } }' >/dev/null 2>&1
  zcat < "$F" | python3 -c "
import sys
from collections import defaultdict
d=defaultdict(dict); hdr=None
for line in sys.stdin:
    if line.startswith('##'): continue
    f=line.split()
    if hdr is None: hdr={k:i for i,k in enumerate(f)}; continue
    if not f[hdr['ID']].startswith('DLGAP1.'): continue
    d[f[hdr['ID']]][f[hdr['TEST']]]=f[hdr['LOG10P']]
for k in sorted(d):
    if 'ADD' in d[k] and 'ADD-SKAT' in d[k]:
        same = d[k]['ADD']==d[k]['ADD-SKAT']
        print('    %-26s burden %-10s SKAT %-10s %s' % (k, d[k]['ADD'], d[k]['ADD-SKAT'], 'IDENTICAL' if same else 'differ'))
"
done

echo
echo "=============== power, MACthr 10 (default) ==============="
python3 scripts/13_power_analyze.py --tag p2_mac10 --design results/power_design_pilot3.tsv \
  --out results/power_mac10.tsv 2>&1 | sed -n '1,40p'

echo
echo "=============== power, MACthr 1 (no collapsing) ==============="
python3 scripts/13_power_analyze.py --tag p2_mac1 --design results/power_design_pilot3.tsv \
  --out results/power_mac1.tsv 2>&1 | sed -n '1,40p'
