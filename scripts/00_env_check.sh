#!/usr/bin/env bash
# Day 0: verify the toolchain before anything downloads.
# Run:  bash scripts/00_env_check.sh
# Paste the whole output back.

set -u

REMOTE_VCF="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz"

echo "=== 1. bcftools ==="
bcftools --version 2>&1 | head -3
echo
echo "--- htslib build features (need libcurl for remote slicing) ---"
bcftools --version 2>&1 | grep -i "Using htslib"
# htsfile prints the compiled-in plugins/features
htsfile --version 2>&1 | head -5 || echo "htsfile not found"

echo
echo "=== 2. plink2 ==="
plink2 --version 2>&1 | head -2 || echo "MISSING: plink2"

echo
echo "=== 3. regenie (need >= 3.2.5 for the --vc-tests syntax used here) ==="
regenie --version 2>&1 | head -3 || echo "MISSING: regenie"

echo
echo "=== 4. tabix / samtools ==="
tabix 2>&1 | head -2
samtools --version 2>&1 | head -2 || echo "samtools not found (optional)"

echo
echo "=== 5. python ==="
python3 --version
python3 -c "import pandas, numpy, scipy; print('pandas', pandas.__version__, '| numpy', numpy.__version__, '| scipy', scipy.__version__)" 2>&1

echo
echo "=== 6. disk space on this volume ==="
df -h . | tail -2

echo
echo "=== 7. THE DECIDING TEST: can bcftools read a remote indexed VCF? ==="
echo "(this streams only the header, a few KB, not the file)"
timeout 120 bcftools view -h "$REMOTE_VCF" 2>&1 | tail -3
RC=$?
echo
if [ $RC -eq 0 ]; then
  echo "RESULT: remote slicing WORKS. We never download a whole chromosome."
else
  echo "RESULT: remote read failed (exit $RC)."
  echo "  If the error mentions 'protocol not supported' or 'unknown scheme',"
  echo "  htslib was built without libcurl. Fix with:  brew reinstall bcftools"
  echo "  If it timed out, it is a network problem, not a build problem."
fi

echo
echo "=== 8. sample count sanity check (expect 3202) ==="
timeout 120 bcftools query -l "$REMOTE_VCF" 2>/dev/null | wc -l
