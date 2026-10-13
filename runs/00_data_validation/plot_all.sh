#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/../.."

for stage in 00_branch_distribution 01_source_comparison 02_gen_distribution 03_gen_source_comparison 04_reco_distribution; do
    bash "runs/00_data_validation/$stage.sh"
done
