#!/usr/bin/env bash
# Day 1, step 3: slice the panel regions out of the 1000G 30x callset.
#
# Per chromosome, resumable. Each chromosome writes data/slice/<chr>.vcf.gz and
# is skipped on rerun if that file already has an index.
#
# Expect roughly 1.5 to 2 GB of download across 17 chromosomes and 4.98 Mb of
# regions, at a few hundred MB per Mb of region with 3202 samples in the source.
#
# Run from the repo root:  bash scripts/03_slice.sh

set -uo pipefail

BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot"
PREFIX="20201028_CCDG_14151_B01_GRM_WGS_2020-08-05"
REF="ref/GRCh38.panel_chroms.fa"
BED="config/panel_regions.bed"
KEEP="config/keep_2504.txt"
OUT="data/slice"
LOG="results/slice_stats.tsv"

mkdir -p "$OUT" results config

echo "bcftools: $(command -v bcftools)  $(bcftools --version 2>&1 | head -1)"
for f in "$REF" "$REF.fai" "$BED" "$KEEP"; do
  [ -s "$f" ] || { echo "MISSING: $f"; exit 1; }
done
echo "regions : $(wc -l < "$BED" | tr -d ' ') intervals, $(awk '{s+=$3-$2} END{printf "%.2f Mb", s/1e6}' "$BED")"
echo "samples : $(wc -l < "$KEEP" | tr -d ' ')"
echo

CHROMS=$(cut -f1 "$BED" | sort -u -V)

[ -s "$LOG" ] || printf "chrom\traw\tpass\tafter_subset\tsplit_biallelic\n" > "$LOG"

for c in $CHROMS; do
  VCF="$OUT/${c}.vcf.gz"
  if [ -s "${VCF}.tbi" ]; then
    echo "$c  cached"
    continue
  fi

  awk -v c="$c" '$1==c' "$BED" > "$OUT/${c}.bed"
  NREG=$(wc -l < "$OUT/${c}.bed" | tr -d ' ')
  URL="${BASE}/${PREFIX}_${c}.recalibrated_variants.vcf.gz"
  echo "$c  $NREG regions ..."
  START=$(date +%s)

  # One pass, in this order:
  #   -R  regions (BED, 0-based half-open)
  #   -S  subset to the 2504 unrelated
  #   -f PASS  keep VQSR-passing sites only
  #   -c 1  drop sites monomorphic AFTER the subset (applied post-subsetting)
  # then split multiallelics and left-align against the callset's own reference,
  # then drop exact duplicates created by the split.
  #
  # No missingness or HWE filter here on purpose. Call rate, HWE departure and
  # heterozygote excess are the measurements that test the mappability
  # prediction in docs/mappability_result.md, so they must be computed on
  # unfiltered calls. Filtering happens at the mask-building step.
  if bcftools view -R "$OUT/${c}.bed" -S "$KEEP" --force-samples \
        -f PASS -c 1 -Ou "$URL" 2>"$OUT/${c}.err" \
     | bcftools norm -m -any -f "$REF" -Ou 2>>"$OUT/${c}.err" \
     | bcftools norm -d exact -Oz -o "$VCF" 2>>"$OUT/${c}.err"
  then
    tabix -f -p vcf "$VCF"
    N=$(bcftools index -n "$VCF" 2>/dev/null)
    SEC=$(( $(date +%s) - START ))
    SZ=$(du -h "$VCF" | cut -f1)
    echo "    $N variants, $SZ, ${SEC}s"
    printf "%s\t\t\t\t%s\n" "$c" "$N" >> "$LOG"
  else
    echo "    FAILED. last lines of $OUT/${c}.err:"
    tail -5 "$OUT/${c}.err"
    rm -f "$VCF"
  fi
done

echo
echo "================== summary =================="
TOT=0
for c in $CHROMS; do
  VCF="$OUT/${c}.vcf.gz"
  if [ -s "${VCF}.tbi" ]; then
    N=$(bcftools index -n "$VCF")
    TOT=$((TOT + N))
    printf "  %-6s %9d variants  %6s\n" "$c" "$N" "$(du -h "$VCF" | cut -f1)"
  else
    printf "  %-6s MISSING\n" "$c"
  fi
done
echo "  ------"
printf "  %-6s %9d variants\n" "total" "$TOT"
echo
echo "  samples in output: $(bcftools query -l "$OUT/$(echo "$CHROMS" | head -1).vcf.gz" 2>/dev/null | wc -l | tr -d ' ')"
echo
echo "DONE. Next: per-variant QC (call rate, HWE, heterozygote excess) to test"
echo "the mappability prediction, then bcftools csq annotation."
