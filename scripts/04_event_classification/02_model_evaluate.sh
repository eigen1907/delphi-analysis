#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Evaluate the saved model and write metrics and diagnostic plots.
CAMPAIGN="${CAMPAIGN:-20261015_jongwon}"
FEATURE_SET="${FEATURE_SET:-detector_pid}"
exec uv run python -B -m event_classification.evaluate \
    --input "output/dataset/$CAMPAIGN" \
    --model "output/models/${CAMPAIGN}_${FEATURE_SET}_default/bdt.joblib" \
    --output "output/plots/$CAMPAIGN/04_event_classification/02_model_evaluate/${FEATURE_SET}_default" \
    "$@"
