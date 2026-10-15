#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

# Output defaults to output/plots/<campaign>/03_sample_validation/<split>.
for plot in branch_distribution source_comparison gen_distribution gen_source_comparison reco_distribution; do
  uv run python -m "sample_validation.overall.$plot" \
    --input output/dataset/20261015_jongwon \
    --split train \
    "$@"
done
