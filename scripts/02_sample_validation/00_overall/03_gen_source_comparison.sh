#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

# Output defaults to output/plots/<campaign>/02_sample_validation/<split>.
# Checks default to output/dataset/checks/<campaign>/<split>/gen_compare/.
exec uv run python -B -m sample_validation.overall.gen_source_comparison \
  --input output/dataset/20260828_florian \
  --split train \
  "$@"
