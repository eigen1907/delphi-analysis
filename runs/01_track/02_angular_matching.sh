#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

dataset=data/20260828_florian
uv run python scripts/01_track/02_angular_matching.py -i "$dataset"
