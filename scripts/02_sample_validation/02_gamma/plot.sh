#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

# Output defaults to output/plots/<campaign>/02_sample_validation/<split>.
exec uv run python -B -m sample_validation.gamma.plot \
  --input output/dataset/20260828_florian \
  --split train \
  --ncpu 5 \
  --stage all \
  "$@"
