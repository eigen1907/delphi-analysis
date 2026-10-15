#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
work=$(cd "$1" && pwd)
nrun=$2
events=$(cat "$work/events.txt")
# A restarted incomplete job is simulated again from its Gen input.
rm -f "$work/.done" "$work/simana.fadsim" "$work/simana.fadana" "$work/simana.sdst"

apptainer exec --cleanenv --no-home --pwd /work \
    --bind "$work:/work" --bind "$CACHE/condition-data:/eos/opendata/delphi/condition-data:ro" \
    "$IMAGE" /bin/bash -c '
    set -e
    exec runsim -VERSION v94c -LABO CERN -NRUN "$1" -EBEAM 45.625 -NEVMAX "$2" -gext my_events.fadgen
    ' bash "$nrun" "$events"
