#!/usr/bin/env bash
set -e

cd "$(dirname "$0")/../.."

uv run python scripts/isr_study/dump_truth_trees.py
