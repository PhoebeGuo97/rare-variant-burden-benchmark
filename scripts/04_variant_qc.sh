#!/usr/bin/env bash
# Day 2, step 4: per-variant QC on the unfiltered slice.
#
# This is the measurement that tests the prediction in docs/mappability_result.md:
# the eight segdup-flagged genes should show heterozygote excess, while MAPT,
# CRHR1 and all of 16p11.2 should look ordinary.
#
# Run:  conda activate rvb && bash scripts/04_variant_qc.sh

set -uo pipefail

SLICE="data/slice"
OUT="results"
mkdir -p "$OUT" config

echo "=================== 0. REGENIE 4.x syntax check ==================="
echo "the plan was written against the 3.x docs, so confirm before relying on it"
regenie --help 2>&1 | grep -i -B1 -A4 "vc-tests" || echo "  no --vc-tests in help output, PASTE THE FULL HELP"
echo
regenie --help 2>&1 | grep -i -E "aaf-bins|build-mask|anno-file|set-list|mask-def|joint|singleton" || true

echo
echo "=================== 1. EUR sample list ==================="
# Heterozygote excess pooled across 5 superpopulations is confounded: real
# population structure causes a heterozygote DEFICIT at differentiated variants
# (Wahlund effect), while collapsed paralogs cause EXCESS. The two push in
# opposite directions, so pooled is still informative, but a within-population
# estimate is the clean one. EUR is used because it is the largest single
# superpopulation at n=503 after AFR, and matches the dissertation's GWAS.
awk -F'\t' '$3=="EUR"{print $1}' config/sample_ancestry.tsv > config/keep_eur.txt
echo "  EUR samples: $(wc -l < config/keep_eur.txt | tr -d ' ')"

echo
echo "=================== 2. pooled QC (all 2504) ==================="
POOLED="$OUT/variant_qc_pooled.tsv"
printf "CHROM\tPOS\tREF\tALT\tAC\tAN\tF_MISSING\tExcHet\tHWE\n" > "$POOLED"
for V in "$SLICE"/chr*.vcf.gz; do
  c=$(basename "$V" .vcf.gz)
  printf "  %-6s " "$c"
  bcftools +fill-tags "$V" -Ou -- -t AN,AC,F_MISSING,ExcHet,HWE 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/F_MISSING\t%INFO/ExcHet\t%INFO/HWE\n' \
      >> "$POOLED" 2>/dev/null
  echo "$(( $(grep -c "^$c	" "$POOLED") )) rows"
done
echo "  total: $(( $(wc -l < "$POOLED") - 1 )) variants"

echo
echo "=================== 3. EUR-only QC ==================="
EURQC="$OUT/variant_qc_eur.tsv"
printf "CHROM\tPOS\tREF\tALT\tAC_EUR\tAN_EUR\tExcHet_EUR\tHWE_EUR\n" > "$EURQC"
for V in "$SLICE"/chr*.vcf.gz; do
  c=$(basename "$V" .vcf.gz)
  printf "  %-6s\n" "$c"
  bcftools view -S config/keep_eur.txt --force-samples -Ou "$V" 2>/dev/null \
    | bcftools +fill-tags -Ou -- -t AN,AC,ExcHet,HWE 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/ExcHet\t%INFO/HWE\n' \
      >> "$EURQC" 2>/dev/null
done
echo "  total: $(( $(wc -l < "$EURQC") - 1 )) variants"

gzip -f "$POOLED" "$EURQC"
echo
echo "wrote ${POOLED}.gz and ${EURQC}.gz"
echo
echo "Next: python3 scripts/04b_qc_summary.py"
