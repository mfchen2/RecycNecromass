#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PICRUST2_BIN="${PICRUST2_BIN:-$PROJECT_DIR/.miniconda/envs/picrust2-osx64/bin}"
THREADS="${THREADS:-8}"

if [[ -x "$PICRUST2_BIN/picrust2_pipeline.py" ]]; then
    export PATH="$PICRUST2_BIN:$PATH"
fi

if ! command -v picrust2_pipeline.py >/dev/null 2>&1; then
    cat >&2 <<EOF
Could not find picrust2_pipeline.py.

Set PICRUST2_BIN to the bin directory of a PICRUSt2 environment, for example:
  PICRUST2_BIN=/path/to/picrust2/bin bash scripts/run_picrust2_ko_predictions.sh
EOF
    exit 1
fi

INPUT_FASTA="${INPUT_FASTA:-$PROJECT_DIR/qiime2-results-local/exports/rep-seqs/dna-sequences.fasta}"
INPUT_BIOM="${INPUT_BIOM:-$PROJECT_DIR/qiime2-results-local/exports/feature-table/feature-table.biom}"
OUTDIR="${OUTDIR:-$PROJECT_DIR/outputs/picrust2_ko_predictions}"

test -f "$INPUT_FASTA"
test -f "$INPUT_BIOM"
mkdir -p "$(dirname "$OUTDIR")"

picrust2_pipeline.py \
    -s "$INPUT_FASTA" \
    -i "$INPUT_BIOM" \
    -o "$OUTDIR" \
    -p "$THREADS" \
    --stratified

Rscript "$PROJECT_DIR/scripts/picrust2_predicted_vs_measured_ko.R"
