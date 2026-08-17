#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RAW_ROOT="${RAW_ROOT:-$PROJECT_DIR/result_X202SC25058062-Z01-F002/01.RawData}"
LOCAL_QIIME_BIN="${LOCAL_QIIME_BIN:-$PROJECT_DIR/.miniconda/envs/qiime2-amplicon-2026.1/bin}"

if [[ -x "$LOCAL_QIIME_BIN/qiime" ]]; then
    export PATH="$LOCAL_QIIME_BIN:$PATH"
    export CONDA_PREFIX="${CONDA_PREFIX:-${LOCAL_QIIME_BIN%/bin}}"
    export R_HOME="${R_HOME:-$CONDA_PREFIX/lib/R}"
fi
mkdir -p "$PROJECT_DIR/.cache/matplotlib"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PROJECT_DIR/.cache/matplotlib}"

if ! command -v qiime >/dev/null 2>&1; then
    cat >&2 <<'EOF'
Could not find the qiime command.

Activate a QIIME 2 environment first, then rerun:
  conda activate qiime2-amplicon-2026.1
  bash scripts/run_qiime2_local.sh

On Apple Silicon, QIIME 2 currently uses the macOS Intel/Rosetta conda build.
EOF
    exit 1
fi

"$PROJECT_DIR/scripts/build_manifest_from_metadata.sh" "$RAW_ROOT"

export PROJECT_DIR
export MANIFEST="${MANIFEST:-$PROJECT_DIR/manifests/manifest-pe.tsv}"
export METADATA="${METADATA:-$PROJECT_DIR/metadata/sample-metadata.tsv}"
export OUTDIR="${OUTDIR:-$PROJECT_DIR/qiime2-results-local}"
export CLASSIFIER="${CLASSIFIER:-$PROJECT_DIR/reference/silva-138-99-nb-classifier.qza}"
if [[ -f "$CLASSIFIER" ]]; then
    export SKIP_TAXONOMY="${SKIP_TAXONOMY:-0}"
else
    export SKIP_TAXONOMY="${SKIP_TAXONOMY:-1}"
fi
export SAMPLING_DEPTH="${SAMPLING_DEPTH:-0}"
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-8}"

bash "$PROJECT_DIR/scripts/run_qiime2_lawrencium.sbatch"
