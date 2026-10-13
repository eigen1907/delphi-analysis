#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/02_gamma/01_sample_distribution.py -i data/20260828_florian
