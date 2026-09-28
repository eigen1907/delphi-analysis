#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/.."

uv run python scripts/plot/isr-fsr-photons.py -i data/20260828_florian
