#!/usr/bin/env bash
# Run `ferro prepare` once to download reference data for ferro-hgvs normalization.
#
# Usage:
#   ./docker/ferro_prepare.sh [output-dir] [extra-ferro-args...]
#
# The output directory is created on the host and mounted into the container.
# Default output dir: ./ferro-reference
#
# Examples:
#   ./docker/ferro_prepare.sh                          # default output dir
#   ./docker/ferro_prepare.sh /data/ferro-reference   # custom output dir
#   ./docker/ferro_prepare.sh ./ferro-ref --genome none --no-refseqgene --no-lrg
#
# The --genome none --no-refseqgene --no-lrg flags produce a minimal download
# (~200MB: transcripts + cdot only) sufficient for coding-variant normalization.
# Omit those flags for the full reference (~6GB) including genomic/intronic support.

set -euo pipefail

OUTPUT_DIR="${1:-./ferro-reference}"
shift || true  # remaining args forwarded to ferro prepare

# Resolve to absolute path before mounting
OUTPUT_DIR="$(cd "$(dirname "$OUTPUT_DIR")" 2>/dev/null && pwd)/$(basename "$OUTPUT_DIR")"
mkdir -p "$OUTPUT_DIR"

IMAGE="weaver-bench"

# Build the image if it doesn't exist
if ! docker image inspect "$IMAGE" > /dev/null 2>&1; then
    echo "Image '$IMAGE' not found — building from docker/Dockerfile.bench..."
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    docker build -f "$SCRIPT_DIR/Dockerfile.bench" -t "$IMAGE" "$SCRIPT_DIR/.."
fi

echo "Running: ferro prepare --output-dir /ferro-reference $*"
echo "Output will be written to: $OUTPUT_DIR"

docker run --rm \
    -v "$OUTPUT_DIR:/ferro-reference" \
    --entrypoint ferro \
    "$IMAGE" \
    prepare --output-dir /ferro-reference "$@"

echo "Done. Reference data written to: $OUTPUT_DIR"
echo ""
echo "To run validation with ferro normalization:"
echo "  docker run --rm \\"
echo "    -v /path/to/data:/data \\"
echo "    -v $OUTPUT_DIR:/ferro-reference \\"
echo "    $IMAGE \\"
echo "    clinvar_baseline_validation.tsv \\"
echo "    --gff GRCh38_latest_genomic.gff.gz \\"
echo "    --fasta GCF_000001405.40_GRCh38.p14_genomic.fna \\"
echo "    --ferro-reference /ferro-reference \\"
echo "    --output-file /data/results.tsv"
