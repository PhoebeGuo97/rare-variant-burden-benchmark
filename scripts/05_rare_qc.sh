#!/usr/bin/env bash
# Day 2, step 5: extract site-level GATK annotations, which are the QC metrics
# that actually apply to rare variants.
#
# Why this exists: heterozygote excess and HWE need expected homozygotes, so
# they have no power below roughly 1/sqrt(N). Every variant entering a mask here
# is MAF < 0.01, far below that. The metrics below are per-site and frequency
# independent, so they reach the variants under test:
#
#   MQ    RMS mapping quality at the site. Segmental duplications drag it down.
#   MQ0   count of mapping-quality-zero reads. Direct read of ambiguous mapping.
#   ME    Mendelian error rate in the 1000G trios. Computed on the full 3202
#         callset before we subset to the 2504 unrelated, so it survives in the
#         INFO field. This is the gold standard error metric and it works at
#         any frequency.
#   VQSLOD, QD, InbreedingCoeff, DP, FS, SOR: the rest of the VQSR annotation set.
#
# These are INFO fields carried from the original callset. bcftools does not
# recompute them on sample subsetting, which is what we want here: ME in
# particular is only meaningful because it was computed while the trios were
# still present.
#
# Run:  conda activate rvb && bash scripts/05_rare_qc.sh

set -uo pipefail

SLICE="data/slice"
OUT="results/site_annot.tsv"
mkdir -p results

echo "bcftools: $(command -v bcftools)"
echo

printf "CHROM\tPOS\tREF\tALT\tMQ\tMQ0\tDP\tQD\tVQSLOD\tInbreedingCoeff\tME\tFS\tSOR\n" > "$OUT"

for V in "$SLICE"/chr*.vcf.gz; do
  c=$(basename "$V" .vcf.gz)
  printf "  %-6s " "$c"
  before=$(wc -l < "$OUT")
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/MQ\t%INFO/MQ0\t%INFO/DP\t%INFO/QD\t%INFO/VQSLOD\t%INFO/InbreedingCoeff\t%INFO/ME\t%INFO/FS\t%INFO/SOR\n' \
    "$V" >> "$OUT" 2>/dev/null
  echo "$(( $(wc -l < "$OUT") - before )) rows"
done

echo
echo "total: $(( $(wc -l < "$OUT") - 1 )) variants"

echo
echo "--- field completeness (fraction that are '.') ---"
awk -F'\t' 'NR>1{n++; for(i=5;i<=13;i++) if($i==".") m[i]++}
     END{split("MQ MQ0 DP QD VQSLOD InbCoeff ME FS SOR",h," ");
         for(i=5;i<=13;i++) printf "  %-12s %6.2f%% missing\n", h[i-4], 100*m[i]/n}' "$OUT"

gzip -f "$OUT"
echo
echo "wrote ${OUT}.gz"
echo
echo "If ME is mostly missing, the Mendelian error route is closed and the"
echo "verdict has to rest on MQ, MQ0 and VQSLOD. Check the line above."
echo
echo "Next: python3 scripts/05b_rare_qc_summary.py"
