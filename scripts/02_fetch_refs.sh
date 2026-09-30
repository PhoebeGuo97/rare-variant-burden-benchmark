#!/usr/bin/env bash
# Day 1, step 2: toolchain report, sample lists, and reference chromosomes.
# Safe to rerun; everything is cached.
#
# Run from the repo root:  bash scripts/02_fetch_refs.sh

set -uo pipefail

FTP="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp"
PANEL_URL="${FTP}/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
TEST_VCF="${FTP}/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz"

mkdir -p ref config data/slice results

echo "================== 0. toolchain =================="
[ -n "${CONDA_DEFAULT_ENV:-}" ] && echo "conda env ACTIVE: $CONDA_DEFAULT_ENV  <-- note which binaries this shadows"
for t in bcftools tabix bgzip plink2 regenie python3; do
  p=$(command -v "$t" 2>/dev/null)
  if [ -z "$p" ]; then
    printf "  %-9s MISSING\n" "$t"
  else
    v=$("$t" --version 2>&1 | head -1)
    printf "  %-9s %-46s %s\n" "$t" "$p" "$v"
  fi
done

echo
echo "  remote VCF read test:"
if bcftools view -h "$TEST_VCF" 2>/dev/null | tail -1 | grep -q "^#CHROM"; then
  NS=$(bcftools query -l "$TEST_VCF" 2>/dev/null | wc -l | tr -d ' ')
  echo "    OK. htslib has libcurl. Samples in file: $NS (expect 3202)"
else
  echo "    FAILED. bcftools cannot read the remote VCF."
  echo "    Most likely htslib was built without libcurl:  brew reinstall bcftools"
fi

echo
echo "================== 1. samples =================="
if [ ! -s config/samples_2504.tsv ]; then
  curl -sSL --fail -o config/samples_2504.tsv "$PANEL_URL" || exit 1
fi
awk 'NR>1 && NF>=3 {print $1}' config/samples_2504.tsv > config/keep_2504.txt
awk 'NR>1 && NF>=3 {print $1"\t"$2"\t"$3}' config/samples_2504.tsv > config/sample_ancestry.tsv
N=$(wc -l < config/keep_2504.txt | tr -d ' ')
echo "  unrelated samples: $N (expect 2504)"
awk '{c[$3]++} END{for(k in c) printf "    %-5s %d\n", k, c[k]}' config/sample_ancestry.tsv | sort

echo
echo "================== 2. reference =================="
# bcftools csq needs the exact reference the callset was built against.
# Fetched by HTTP range request from the remote .fai, so no samtools needed.
CHROMS=$( (cut -f1 config/panel_regions.bed; echo chr20; echo chr21; echo chr22) | sort -u -V | tr '\n' ' ')
echo "  chromosomes: $CHROMS"
REF="ref/GRCh38.panel_chroms.fa"

if [ -s "${REF}.fai" ]; then
  echo "  using cached $REF"
else
  python3 scripts/02_fetch_ref.py --chroms $CHROMS --out "$REF" || exit 1
fi

echo
awk '{s+=$2} END{printf "  %d contigs, %.2f Gb\n", NR, s/1e9}' "${REF}.fai"

echo
echo "DONE."
