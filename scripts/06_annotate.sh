#!/usr/bin/env bash
# Day 2, step 6: build the panel GFF3, run bcftools csq, set variant IDs.
#
# Run:  conda activate rvb && bash scripts/06_annotate.sh

set -uo pipefail

REL=${REL:-112}
GFF_SRC="ref/Homo_sapiens.GRCh38.${REL}.chr.gff3.gz"
GFF="ref/panel.chr.gff3.gz"
REF="ref/GRCh38.panel_chroms.fa"
SLICE="data/slice"
OUT="data/csq"
mkdir -p "$OUT" results

echo "=============== 0. REGENIE 4.x option names ==============="
echo "confirming before the mask files are written to a 3.x-era format"
regenie --help 2>&1 | grep -i -E "anno-file|set-list|mask-def|aaf-bins|aaf-file|build-mask|vc-tests|singleton|write-mask" || \
  echo "  nothing matched, PASTE the full 'regenie --help' output"

echo
echo "=============== 1. panel GFF3 ==============="
if [ -s "$GFF" ]; then
  echo "  using cached $GFF"
else
  python3 scripts/06_make_gff.py --gff "$GFF_SRC" --out ref/panel.chr.gff3 || exit 1
  bgzip -f ref/panel.chr.gff3
  echo "  size: $(du -h "$GFF" | cut -f1)"
fi

echo
echo "=============== 2. bcftools csq ==============="
# --local-csq: consequences per variant, not haplotype aware. The callset is
# unphased, and masks are built per variant anyway, so compound-heterozygote
# reasoning would add nothing and would require phase we do not have.
# --set-id: the raw callset has "." for nearly every ID. REGENIE matches
# variants between the annotation file, the set list and the genotypes by ID,
# so they have to be set once, here, and never regenerated differently later.
for V in "$SLICE"/chr*.vcf.gz; do
  c=$(basename "$V" .vcf.gz)
  O="$OUT/${c}.csq.vcf.gz"
  if [ -s "${O}.tbi" ]; then
    echo "  $c  cached"
    continue
  fi
  printf "  %-6s " "$c"
  if bcftools csq -f "$REF" -g "$GFF" --local-csq --ncsq 30 -Ou "$V" 2>"$OUT/${c}.err" \
     | bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Oz -o "$O" 2>>"$OUT/${c}.err"
  then
    tabix -f -p vcf "$O"
    N=$(bcftools index -n "$O")
    NB=$(bcftools query -f '%INFO/BCSQ\n' "$O" 2>/dev/null | grep -vc '^\.$')
    echo "$N variants, $NB with a consequence"
  else
    echo "FAILED"
    tail -5 "$OUT/${c}.err"
  fi
done

echo
echo "=============== 3. consequence spectrum ==============="
# What csq actually produced, before any of it is trusted.
for V in "$OUT"/chr*.csq.vcf.gz; do
  bcftools query -f '%INFO/BCSQ\n' "$V" 2>/dev/null
done \
| tr ',' '\n' | sed 's/^[@*]//' | cut -d'|' -f1 | tr '&' '\n' \
| grep -v '^\.$' | grep -v '^$' | sort | uniq -c | sort -rn \
| awk '{printf "  %-28s %d\n", $2, $1}'

echo
echo "DONE. Next: python3 scripts/06b_make_masks.py"
