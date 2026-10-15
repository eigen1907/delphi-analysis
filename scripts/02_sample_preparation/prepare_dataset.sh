#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Pair native/raw events and create train, validation and test files per sample.
exec uv run python -m sample_preparation.dataset \
    --input data/20261015_jongwon \
    --counts 60000 15000 15000 \
    --seed 1907 \
    "$@"
