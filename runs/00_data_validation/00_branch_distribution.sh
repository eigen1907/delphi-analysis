#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/../.."

uv run python scripts/00_data_validation/00_branch_distribution.py \
    -i data/202606xx_jongwon/dataset/20260606_100kTest
