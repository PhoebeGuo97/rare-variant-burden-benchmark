#!/usr/bin/env bash
# Day 2/3, step 10: the full simulation grid.  v2.
#
# v3: adds an ARMS selector so a partly-successful v1 run can be topped up.
# v1 died at phenotype 246 with "cannot write file". REGENIE opens one output
# file per phenotype and holds them all open at once, and the macOS default
# open-file limit is 256. Two fixes, belt and braces:
#   1. raise ulimit -n at the top of this script
#   2. run phenotypes in chunks, so the number of simultaneously open files is
#      bounded by CHUNK no matter what the limit happens to be
#
# Chunking also makes resume finer grained: a crash costs one chunk, not an arm.
# It costs re-reading the genotypes once per chunk, which is a minute or so.
#
# Run:  conda activate rvb && nohup bash scripts/10_run_full_v2.sh > results/full_run.out 2>&1 &
#       tail -f results/full_run.out

set -uo pipefail

# ---------------------------------------------------------------- file limits
echo "open file limit before: $(ulimit -n)"
ulimit -n 10240 2>/dev/null || ulimit -n 4096 2>/dev/null || ulimit -n 2048 2>/dev/null || true
echo "open file limit after : $(ulimit -n)"
echo "hard limit            : $(ulimit -Hn)"
echo

P="data/plink/panel"
R="results/regenie"
COV="config/covariates.tsv"
PCS="PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
mkdir -p "$R" results

N_QT=${N_QT:-500}
N_BT=${N_BT:-200}
N_ST=${N_ST:-100}
CHUNK=${CHUNK:-100}

# Which arms to run. Default all. Set ARMS to fill in only the ones that failed,
# e.g. ARMS="a1 b1 b2" after the v1 run completed a2 through a5 successfully.
ARMS=${ARMS:-"a1 a2 a3 a4 a5 b1 b2"}
want() { [[ " $ARMS " == *" $1 "* ]]; }
echo "arms to run: $ARMS"

COMMON=(
  --step 2
  --pgen "$P"
  --anno-file config/regenie_annot.tsv
  --set-list config/regenie_setlist.tsv
  --mask-def config/regenie_masks.def
  --aaf-bins 0.001,0.005,0.01
  --vc-tests skat,skato,skato-acat,acatv,acato-full
  --build-mask max
  --minMAC 1
  --ignore-pred
  --bsize 200
  --gz
)

# run_arm <tag> <phenoFile> [extra regenie args...]
# Splits the phenotype file's columns into chunks of CHUNK and runs each.
run_arm() {
  local tag=$1; shift
  local pheno=$1; shift
  if [ ! -s "$pheno" ]; then
    echo "  MISSING phenotype file: $pheno"
    return 1
  fi
  # columns 3.. are the traits; columns 1-2 are FID and IID
  local -a cols
  IFS=$'\t' read -r -a hdr < "$pheno"
  cols=("${hdr[@]:2}")
  local n=${#cols[@]}
  local nchunk=$(( (n + CHUNK - 1) / CHUNK ))
  echo "[$(date +%H:%M:%S)] $tag : $n traits in $nchunk chunk(s) of $CHUNK"

  local i=0 c=0
  while [ $i -lt $n ]; do
    c=$((c + 1))
    local sub=("${cols[@]:$i:$CHUNK}")
    local list
    list=$(IFS=,; echo "${sub[*]}")
    local sub_tag="${tag}_c${c}"
    if [ -s "$R/${sub_tag}.log" ] && grep -q "End time" "$R/${sub_tag}.log" 2>/dev/null; then
      echo "  [$(date +%H:%M:%S)] $sub_tag  cached"
      i=$((i + CHUNK)); continue
    fi
    local t0=$(date +%s)
    regenie "${COMMON[@]}" --phenoFile "$pheno" --phenoColList "$list" "$@" \
      --out "$R/$sub_tag" > "$R/${sub_tag}.stdout" 2>&1
    local rc=$?
    local nf
    nf=$(ls "$R/${sub_tag}"_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ')
    echo "  [$(date +%H:%M:%S)] $sub_tag  exit $rc  $(( $(date +%s) - t0 ))s  ${nf}/${#sub[@]} files"
    if [ $rc -ne 0 ] || [ "$nf" -ne "${#sub[@]}" ]; then
      echo "    INCOMPLETE, last lines:"
      tail -8 "$R/${sub_tag}.stdout" | sed 's/^/      /'
      rm -f "$R/${sub_tag}.log"
    fi
    i=$((i + CHUNK))
  done
}

echo "=============== phenotypes ==============="
python3 scripts/09_make_phenos.py --n-rep "$N_QT" --tag qt || exit 1
python3 scripts/09_make_phenos.py --n-rep "$N_BT" --tag bt || exit 1
python3 scripts/09_make_phenos.py --n-rep "$N_ST" --tag st || exit 1

echo
echo "=============== ARM A: unstructured null ==============="
want a1 && run_arm a1_qt config/pheno_qt_null_qt.tsv --covarFile "$COV" --covarColList "$PCS"

# Uncorrected score tests on unbalanced binary traits. Expected to fail, and the
# failure is the deliverable, not a bug to fix.
want a2 && run_arm a2_bt10_nofirth config/pheno_bt10_null_bt.tsv --bt \
    --covarFile "$COV" --covarColList "$PCS"
want a3 && run_arm a3_bt01_nofirth config/pheno_bt01_null_bt.tsv --bt \
    --covarFile "$COV" --covarColList "$PCS"

# Same traits, Firth corrected. The paired comparison is the result.
want a4 && run_arm a4_bt10_firth config/pheno_bt10_null_bt.tsv --bt --firth --approx --pThresh 0.01 \
    --covarFile "$COV" --covarColList "$PCS"
want a5 && run_arm a5_bt01_firth config/pheno_bt01_null_bt.tsv --bt --firth --approx --pThresh 0.01 \
    --covarFile "$COV" --covarColList "$PCS"

echo
echo "=============== ARM B: structure confounded null ==============="
want b1 && run_arm b1_struct_nopc config/pheno_qt_struct_st.tsv
want b2 && run_arm b2_struct_pc   config/pheno_qt_struct_st.tsv --covarFile "$COV" --covarColList "$PCS"

echo
echo "=============== inventory ==============="
for tag in a1_qt a2_bt10_nofirth a3_bt01_nofirth a4_bt10_firth a5_bt01_firth \
           b1_struct_nopc b2_struct_pc; do
  n=$(ls "$R/${tag}"_c*_*.regenie.gz 2>/dev/null | wc -l | tr -d ' ')
  printf "  %-20s %5s result files\n" "$tag" "$n"
done

echo
echo "DONE."
