#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Set CAMPAIGN=20260828_florian to use the other dataset.
CAMPAIGN="${CAMPAIGN:-20261015_jongwon}"
FEATURE_SET="${FEATURE_SET:-detector_pid}"
# Both aligned trees provide features; val.root controls early stopping.
exec uv run python -m event_classification.train \
    --input "output/dataset/$CAMPAIGN" \
    --output "output/models/${CAMPAIGN}_${FEATURE_SET}_default" \
    --feature-config "config/pid/features/$FEATURE_SET.json" \
    --hyperparameters config/pid/hyperparameters/default.json \
    --seed 1907 \
    --ncpu 32 \
    "$@"
