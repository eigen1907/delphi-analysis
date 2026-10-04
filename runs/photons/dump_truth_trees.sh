#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/photons/dump-truth-trees.py
