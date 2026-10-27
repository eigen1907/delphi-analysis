#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Apply to each class's test split, preserving the event keys in predictions.
CAMPAIGN="${CAMPAIGN:-20261015_jongwon}"
FEATURE_SET="${FEATURE_SET:-detector_pid}"
exec uv run python -B -m event_classification.apply \
    --input "output/dataset/$CAMPAIGN" \
    --model "output/models/${CAMPAIGN}_${FEATURE_SET}_default/bdt.joblib" \
    --output "output/dataset/$CAMPAIGN/predictions/${FEATURE_SET}_default" \
    --split test \
    "$@"
