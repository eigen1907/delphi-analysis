#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

sample=20260606_100kTest
data_dir=data/202606xx_jongwon
model="plots/bdt/${sample}_pid_standard/bdt.joblib"

for class in Zee Zmumu ZKK Zpipi; do
    uv run python scripts/bdt/apply.py \
        -i "$data_dir/dataset/$sample/$class/nanoaod_raw_sdst.root" \
        -m "$model" \
        -o "$data_dir/ml/${sample}_pid_standard/prediction_${class}.root"
done
