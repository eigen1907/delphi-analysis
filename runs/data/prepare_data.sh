#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

sample=20260606_100kTest
data_dir=data/202606xx_jongwon

uv run python scripts/data/prepare-chunks.py \
    -i "$data_dir/raw/$sample" \
    -o "$data_dir/chunk/$sample"

uv run python scripts/data/build-dataset.py \
    -i "$data_dir/chunk/$sample" \
    -o "$data_dir/dataset/$sample"
