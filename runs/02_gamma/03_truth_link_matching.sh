#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/02_gamma/03_truth_link_matching.py -i data/20260828_florian
