#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/../.."

dataset=data/20260828_florian
uv run python scripts/03_rich/02_angular_matching.py -i "$dataset" -o plots
