#!/usr/bin/env bash
# Day 3, step 15: install Nextflow and make a tiny phenotype file for a real
# end-to-end test of the pipeline.
#
# No meaningful CPU: safe to run while the power grid is still going. Do NOT
# run the pipeline itself until the grid finishes.
#
# Run:  conda activate rvb && bash scripts/15_setup_nextflow.sh

set -uo pipefail

echo "=============== java ==============="
if command -v java >/dev/null 2>&1; then
  java -version 2>&1 | head -2 | sed 's/^/  /'
else
  echo "  MISSING. Nextflow needs Java 17 or newer."
  echo "  brew install --cask temurin"
  exit 1
fi

echo
echo "=============== nextflow ==============="
mkdir -p "$HOME/bin"
if command -v nextflow >/dev/null 2>&1; then
  echo "  already installed: $(command -v nextflow)"
elif [ -x "$HOME/bin/nextflow" ]; then
  echo "  found at \$HOME/bin/nextflow"
else
  echo "  installing to \$HOME/bin ..."
  ( cd "$HOME/bin" && curl -s https://get.nextflow.io | bash ) || exit 1
fi
export PATH="$HOME/bin:$PATH"
nextflow -version 2>&1 | sed 's/^/  /'

echo
echo "  add this to ~/.zshrc if it is not already there:"
echo "    export PATH=\"\$HOME/bin:\$PATH\""

echo
echo "=============== tiny phenotype file ==============="
# 8 traits, drawn from the existing null file so nothing new is simulated.
# Two chunks at chunk_size 4, which is the point: it exercises the chunking and
# the aggregation, not just a single REGENIE call.
python3 - <<'PY'
src = "config/pheno_qt_null_qt.tsv"
dst = "config/pheno_smoke.tsv"
with open(src) as fh, open(dst, "w") as out:
    hdr = fh.readline().rstrip("\n").split("\t")
    keep = [0, 1] + list(range(2, 10))
    out.write("\t".join(hdr[i] for i in keep) + "\n")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        out.write("\t".join(f[i] for i in keep) + "\n")
print("  wrote %s with %d traits" % (dst, 8))
PY

echo
echo "=============== ready ==============="
cat <<'EOF'
  Once the power grid has finished, run:

    cd nextflow
    nextflow run main.nf -entry assoc_only -profile local \
        --pheno ../config/pheno_smoke.tsv \
        --chunk_size 4 \
        --outdir ../results_nf_smoke

  Expect two ASSOC tasks (8 traits / chunk_size 4) and one CALIBRATION task.
  A few minutes. Success means the chunking, the REGENIE invocation and the
  aggregation all work under Nextflow, which is the part worth claiming.

  The -profile local run uses the rvb conda binaries and no containers, so it
  does NOT test the container tags in nextflow.config. Those remain unverified.
EOF
