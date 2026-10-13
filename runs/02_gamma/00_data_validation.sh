#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/02_gamma/00_data_validation.py -i data/20260828_florian
