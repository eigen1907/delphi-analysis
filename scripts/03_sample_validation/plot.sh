#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."

# Output defaults to output/plots/<campaign>/03_sample_validation/<split>.
for study in 00_overall 01_track 02_gamma 03_rich; do
  bash "scripts/03_sample_validation/$study/plot.sh" "$@"
done
