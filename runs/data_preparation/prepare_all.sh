#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/../.."

bash runs/data_preparation/00_prepare_chunks.sh
bash runs/data_preparation/01_build_dataset.sh
