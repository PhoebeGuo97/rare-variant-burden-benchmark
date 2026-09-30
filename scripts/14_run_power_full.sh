#!/usr/bin/env bash
# Day 3, step 14: the full power grid.
#
# 4 genes spanning 17 to 204 causal variants, 3 effect sizes, 2 direction arms,
# 50 replicates = 1200 traits, run twice (MACthr 10 and MACthr 1).
#
# The pilot established the finding on one gene at 25 reps. This exists to show
# it holds across variant counts and to tighten the estimates: 50 reps gives a
# standard error near 0.07 per cell instead of 0.10.
#
# Long. Chunked and resumable, same as the null grid.
#
# Run:  conda activate rvb && nohup bash scripts/14_run_power_full.sh > results/power_full.out 2>&1 &
#       tail -20 results/power_full.out        # in a DIFFERENT tab, and not with -f

set -uo pipefail
echo "open file limit: $(ulimit -n) -> "
ulimit -n 10240 2>/dev/null || true
echo "$(ulimit -n)"

P="data/plink/panel"
R="results/regenie"
COV="config/covariates.tsv"
PCS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
CHUNK=${CHUNK:-100}
mkdir -p "$R" results

echo
echo "=============== simulate ==============="
python3 scripts/12_sim_power.py --mode full || exit 1
PHENO=config/pheno_power_full.tsv

run_arm() {   # run_arm <tag> <macthr>
  local tag=$1 mac=$2
  local -a cols
  IFS=$'\t' read -r -a hdr < "$PHENO"
  cols=("${hdr[@]:2}")
  local n=${#cols[@]}
  echo
  echo "[$(date +%H:%M:%S)] $tag : $n traits, MACthr $mac, chunks of $CHUNK"
  local i=0 c=0
  while [ $i -lt $n ]; do
    c=$((c + 1))
    local sub=("${cols[@]:$i:$CHUNK}")
    local list; list=$(IFS=,; echo "${sub[*]}")
    local st="${tag}_c${c}"
    if [ -s "$R/${st}.log" ] && grep -q "End time" "$R/${st}.log" 2>/dev/null; then
      echo "  [$(date +%H:%M:%S)] $st cached"; i=$((i + CHUNK)); continue
    fi
    local t0=$(date +%s)
    regenie --step 2 --pgen "$P" \
      --phenoFile "$PHENO" --phenoColList "$list" \
      --covarFile "$COV" --covarColList "$PCS" \
      --anno-file config/regenie_annot.tsv \
      --set-list config/regenie_setlist.tsv \
      --mask-def config/regenie_masks.def \
      --aaf-bins 0.001,0.005,0.01 \
      --vc-tests skat,skato,skato-acat,acatv,acato-full \
      --vc-maxAAF 0.01 --vc-MACthr "$mac" \
      --build-mask max --minMAC 1 --ignore-pred --bsize 200 --gz \
      --out "$R/$st" > "$R/${st}.stdout" 2>&1
    local rc=$? nf
    nf=$(ls "$R/${st}"_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ')
    echo "  [$(date +%H:%M:%S)] $st exit $rc $(( $(date +%s) - t0 ))s ${nf}/${#sub[@]} files"
    if [ $rc -ne 0 ] || [ "$nf" -ne "${#sub[@]}" ]; then
      tail -6 "$R/${st}.stdout" | sed 's/^/      /'; rm -f "$R/${st}.log"
    fi
    i=$((i + CHUNK))
  done
}

run_arm f_mac10 10
run_arm f_mac1 1

echo
echo "=============== analyze ==============="
python3 scripts/13_power_analyze.py --tag f_mac10 \
  --design results/power_design_full.tsv --out results/power_full_mac10.tsv
python3 scripts/13_power_analyze.py --tag f_mac1 \
  --design results/power_design_full.tsv --out results/power_full_mac1.tsv

echo
echo "DONE."
