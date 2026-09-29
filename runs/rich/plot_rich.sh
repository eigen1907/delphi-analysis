#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

dataset=data/202606xx_jongwon/dataset/20260606_100kTest

uv run python scripts/rich/rich.py -i "$dataset"
