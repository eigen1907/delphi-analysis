#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/02_gamma/plot.py -i data/20260828_florian
