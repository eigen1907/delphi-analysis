#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../../.."

# Output defaults to output/plots/<campaign>/03_sample_validation/<split>.
exec uv run python -m sample_validation.track.plot \
  --input output/dataset/20261015_jongwon \
  --split train \
  --stage all \
  "$@"
