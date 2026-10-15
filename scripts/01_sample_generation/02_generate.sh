#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
work=$(cd "$1" && pwd)
sample=$2
seed=$3
events=$(cat "$work/events.txt")
gen_events=$(( (events * 10 + 8) / 9 ))
rm -f "$work/.done"
final="$(dirname "$(dirname "$work")")/final_root/$(basename "$work")"
mkdir -p "$final"
cp "$SAMPLE/config/$sample.cmnd" "$work/generator.cmnd"

# Extra Gen records keep the external generator reader away from EOF.
echo "PYTHIA seed: $seed"
apptainer exec --cleanenv --no-home --pwd /work --bind "$work:/work" "$IMAGE" /bin/bash -c '
    set -e
    pythia8_generate "$1" /work/generator.cmnd "$2"
    mv fort.26 my_events.fadgen
    ' bash "$gen_events" "$seed"
mv "$work/pythiastandalone.root" "$final/pythiastandalone.root"
