#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

sample=20260606_100kTest
data_dir=data/202606xx_jongwon

uv run python scripts/bdt/prepare.py \
    -i "$data_dir/chunk/$sample" \
    -o "$data_dir/ml/${sample}_pid"

uv run python scripts/bdt/train.py \
    -i "$data_dir/ml/${sample}_pid" \
    -o "plots/bdt/${sample}_pid_standard"
