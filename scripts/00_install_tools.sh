#!/usr/bin/env bash
# Install the toolchain into a dedicated conda env, not (base).
#
# (base) holds hail with numpy<2, pandas<1.5, scipy<1.10 pins. Installing
# genomics tools there risks a solve that moves those. A separate env also makes
# "which binary am I running" unambiguous, which matters because the toolchain
# report showed tools resolving from several places.
#
# Run:  bash scripts/00_install_tools.sh
# Then: conda activate rvb

set -uo pipefail

ENV_NAME=${ENV_NAME:-rvb}
PKGS="bcftools samtools plink2 regenie"

if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then
  SOLVER=conda
else
  echo "No conda or mamba on PATH. Stop here and say so."
  exit 1
fi
echo "using $SOLVER"
echo "target env: $ENV_NAME"
echo "packages: $PKGS"
echo

# ---- attempt 1: native arm64 ------------------------------------------------
echo "=== attempt 1: native osx-arm64 ==="
$SOLVER create -n "$ENV_NAME" -c conda-forge -c bioconda -y $PKGS
RC=$?

# ---- attempt 2: osx-64 under Rosetta ---------------------------------------
# bioconda does not build every tool for Apple Silicon. regenie is the likely
# gap. An osx-64 env runs fine under Rosetta 2 and is a normal thing to do.
if [ $RC -ne 0 ]; then
  echo
  echo "=== attempt 1 failed. attempt 2: osx-64 under Rosetta ==="
  $SOLVER env remove -n "$ENV_NAME" -y >/dev/null 2>&1
  CONDA_SUBDIR=osx-64 $SOLVER create -n "$ENV_NAME" -c conda-forge -c bioconda -y $PKGS
  RC=$?
  if [ $RC -eq 0 ]; then
    conda config --env --set subdir osx-64 2>/dev/null
    echo "env pinned to osx-64 so later installs stay consistent"
  fi
fi

if [ $RC -ne 0 ]; then
  echo
  echo "Both attempts failed. Paste the solver error."
  echo "Fallback is official binaries:"
  echo "  plink2  https://www.cog-genomics.org/plink/2.0/   (macOS arm64 build)"
  echo "  regenie https://github.com/rgcgithub/regenie/releases"
  exit 1
fi

# ---- verify -----------------------------------------------------------------
echo
echo "=== versions in $ENV_NAME ==="
CONDA_BASE=$(conda info --base)
# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

for t in bcftools samtools plink2 regenie; do
  p=$(command -v "$t" 2>/dev/null)
  if [ -z "$p" ]; then
    printf "  %-9s MISSING\n" "$t"
  else
    printf "  %-9s %s\n" "$t" "$("$t" --version 2>&1 | head -1)"
  fi
done

echo
echo "=== regenie version gate ==="
RV=$(regenie --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "  detected: ${RV:-unknown}"
echo "  need >= 3.2.5 for the --vc-tests syntax this project uses"
if [ -n "$RV" ]; then
  python3 - "$RV" <<'PY'
import sys
v = tuple(int(x) for x in sys.argv[1].split("."))
print("  OK" if v >= (3, 2, 5) else "  TOO OLD, the SKAT-O/ACAT calls will need rewriting")
PY
fi

echo
echo "=== remote read test with this env's bcftools ==="
U="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz"
if bcftools view -h "$U" 2>/dev/null | tail -1 | grep -q "^#CHROM"; then
  echo "  OK, htslib here has libcurl"
else
  echo "  FAILED. This env's bcftools cannot read remote files."
  echo "  Keep using the system bcftools for slicing if that one worked."
fi

echo
echo "DONE. Use this env for every later step:  conda activate $ENV_NAME"
