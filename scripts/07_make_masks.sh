#!/usr/bin/env bash
# Day 2, step 7: extract the csq annotations and build the REGENIE mask files.
#
# The extraction is done by bcftools rather than in python: the annotated VCFs
# are 2.3 GB compressed and carry 2504 genotypes per line, so parsing them in
# python would read tens of GB to reach four fields.
#
# Run:  conda activate rvb && bash scripts/07_make_masks.sh

set -uo pipefail

CSQ="data/csq"
TAB="results/csq_table.tsv"
mkdir -p results config

if [ -s "${TAB}.gz" ]; then
  echo "using cached ${TAB}.gz"
else
  printf "ID\tCHROM\tPOS\tBCSQ\n" > "$TAB"
  for V in "$CSQ"/chr*.csq.vcf.gz; do
    c=$(basename "$V" .csq.vcf.gz)
    printf "  %-6s " "$c"
    before=$(wc -l < "$TAB")
    bcftools query -f '%ID\t%CHROM\t%POS\t%INFO/BCSQ\n' "$V" >> "$TAB" 2>/dev/null
    echo "$(( $(wc -l < "$TAB") - before )) rows"
  done
  gzip -f "$TAB"
fi

echo
echo "=== truncation check ==="
# bcftools warned about too many consequences at chr17:45991586. That warning is
# about the per-sample FORMAT/BCSQ bitmask, which we do not use. This checks
# whether INFO/BCSQ itself was capped, which WOULD matter: a gene's consequence
# could be dropped from a crowded site and silently leave the mask.
python3 - <<'PY'
import gzip
from collections import Counter
mx=0; site=None; hist=Counter()
with gzip.open('results/csq_table.tsv.gz','rt') as fh:
    fh.readline()
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if len(f)<4 or f[3] in ('.',''): continue
        n=f[3].count(',')+1
        hist[n]+=1
        if n>mx: mx=n; site="%s:%s"%(f[1],f[2])
print("  max INFO/BCSQ entries at one site: %d  (at %s)" % (mx, site))
top=sorted(hist)[-5:]
print("  sites at the top counts:", {k:hist[k] for k in top})
print("  --ncsq was 30, giving 62 per-sample slots. If the max above sits")
print("  exactly at a round cap and many sites pile up there, INFO was capped")
print("  too and csq needs rerunning with a higher --ncsq.")
PY

echo
echo "=== building mask files ==="
python3 scripts/07b_make_masks.py
