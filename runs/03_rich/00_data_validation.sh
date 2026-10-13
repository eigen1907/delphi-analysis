#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/../.."

dataset=data/20260828_florian
uv run python scripts/03_rich/00_data_validation.py -i "$dataset" -o plots
