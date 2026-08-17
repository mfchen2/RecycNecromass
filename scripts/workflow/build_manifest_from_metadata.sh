#!/bin/bash
set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
METADATA="${METADATA:-$PROJECT_DIR/metadata/sample-metadata.tsv}"
MANIFEST="${MANIFEST:-$PROJECT_DIR/manifests/manifest-pe.tsv}"
RAW_ROOT="${1:-${RAW_ROOT:-$PROJECT_DIR/result_X202SC25058062-Z01-F002/01.RawData}}"

if [[ "$RAW_ROOT" != /* ]]; then
    echo "RAW_ROOT must be an absolute path: $RAW_ROOT" >&2
    exit 1
fi

mkdir -p "$(dirname "$MANIFEST")"

{
    printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n'
    awk -v raw_root="$RAW_ROOT" 'BEGIN { OFS = "\t" }
        NR > 2 {
            sub(/\r$/, "")
            sample_id = $1
            print sample_id, raw_root "/" sample_id "/" sample_id "_1.fastq.gz", raw_root "/" sample_id "/" sample_id "_2.fastq.gz"
        }' "$METADATA"
} > "$MANIFEST"

echo "Wrote $MANIFEST using RAW_ROOT=$RAW_ROOT"

