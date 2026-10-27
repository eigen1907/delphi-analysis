#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

# Output defaults to output/plots/<campaign>/02_sample_validation/<split>.
for plot in branch_distribution source_comparison gen_distribution gen_source_comparison reco_distribution; do
  uv run python -B -m "sample_validation.overall.$plot" \
    --input output/dataset/20260828_florian \
    --split train \
    "$@"
done
