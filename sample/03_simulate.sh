#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
work=$(cd "$1" && pwd)
nrun=$2
events=$(cat "$work/events.txt")
# A restarted incomplete job is simulated again from its Gen input.
rm -f "$work/.done" "$work/simana.fadsim" "$work/simana.fadana" "$work/simana.sdst"

docker run --init --rm --platform linux/arm64 --entrypoint /bin/bash -w /work \
    -v "$work:/work" -v "$CACHE/condition-data:/eos/opendata/delphi/condition-data:ro" \
    "$SIM_IMAGE" -c '
    set -e
    export DELPHI_INSTALL_DIR=/delphi
    source /delphi/setup.sh
    export DELPHI_DDB=/eos/opendata/delphi/condition-data
    exec runsim -VERSION v94c -LABO CERN -NRUN "$1" -EBEAM 45.625 -NEVMAX "$2" -gext my_events.fadgen
    ' bash "$nrun" "$events"
