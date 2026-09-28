#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_root"

if [[ $# -lt 1 || $# -gt 2 ]]; then
    printf 'Usage: %s SAMPLE_SET [CHUNK_ROOT]\n' "$0" >&2
    exit 2
fi

sample_set="$1"
chunk_root="${2:-data/chunk/$sample_set}"
dataset_root="data/ml/${sample_set}_pid"
model_root="data/models/${sample_set}_pid_standard"
plot_root="plots/bdt/${sample_set}_pid_standard"
log_path="logs/bdt/${sample_set}_pid_standard.log"

mkdir -p "$(dirname "$log_path")"

{
    uv run --locked python scripts/bdt/prepare.py \
        -i "$chunk_root" \
        -o "$dataset_root" \
        --feature-config config/bdt/features/pid.json

    uv run --locked python scripts/bdt/train.py \
        -i "$dataset_root" \
        -o "$model_root" \
        --plot-output "$plot_root" \
        --hyperparameters config/bdt/hyperparameters/standard.json
} 2>&1 | tee "$log_path"
