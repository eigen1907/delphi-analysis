#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

dataset=data/202606xx_jongwon/dataset/20260606_100kTest

for script in branches-all branches-compare gen-check gen-compare reco-check rich gen-reco-track-match-cut gen-reco-track-match-result; do
    uv run python "scripts/plot/$script.py" -i "$dataset"
done
