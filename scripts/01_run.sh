#!/usr/bin/env bash
# Day 1, step 1: panel regions and measured mappability.
# Run from the repo root:  bash scripts/01_run.sh

set -euo pipefail

REL=${REL:-112}
GFF="ref/Homo_sapiens.GRCh38.${REL}.chr.gff3.gz"
GFF_URL="https://ftp.ensembl.org/pub/release-${REL}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${REL}.chr.gff3.gz"
SEGDUP="ref/genomicSuperDups.txt.gz"
SEGDUP_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz"

mkdir -p ref config docs results

if [ ! -s "$GFF" ]; then
  echo "downloading Ensembl ${REL} GFF3 (~50 MB)..."
  curl -L --fail -o "$GFF" "$GFF_URL"
else
  echo "using cached $GFF"
fi

if [ ! -s "$SEGDUP" ]; then
  echo "downloading UCSC hg38 segmental duplications (~10 MB)..."
  curl -L --fail -o "$SEGDUP" "$SEGDUP_URL"
else
  echo "using cached $SEGDUP"
fi

echo
echo "================ 1. panel regions ================"
python3 scripts/01_make_regions.py --gff "$GFF"

echo
echo "================ 1b. measured mappability ================"
python3 scripts/01b_mappability.py

echo
echo "--- slicing intervals (first 10 of $(wc -l < config/panel_regions.bed)) ---"
head -10 config/panel_regions.bed
