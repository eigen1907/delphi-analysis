#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Pair native/raw events and create train, validation and test files per sample.
exec uv run python -B -m sample_preparation.dataset \
    --input data/20260828_florian \
    --seed 1907 \
    "$@"
