#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/isr_study/plot.py -i data/20260828_florian
