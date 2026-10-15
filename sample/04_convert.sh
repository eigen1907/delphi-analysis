#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"
work=$(cd "$1" && pwd)
cd "$PROJECT"
final="$(dirname "$(dirname "$work")")/final_root/$(basename "$work")"
events=$(cat "$work/events.txt")
rm -f "$work/.done"
mkdir -p "$final"
cp "$NANOAOD/config/delphi-nanoaod.yaml" "$work/nanoaod.yaml"

for source in sdst fadana; do
    native=nanoaod.root
    if [[ "$source" == fadana ]]; then native=nanoaod_fadana.root; fi
    mkdir -p "$work/$source"
    printf 'FILE=/input/simana.%s\n' "$source" > "$work/$source/PDLINPUT_SOURCE"
    docker run --init --rm --platform linux/amd64 \
        -v "$work:/input:ro" -v "$work/$source:/work" -v "$final:/out" \
        -v "$CACHE/build:/build:ro" -v "$CACHE/condition-data:/eos/opendata/delphi/condition-data:ro" \
        "$IMAGE" -c '
        set -e
        /build/delphi-nanoaod/delphi-nanoaod -P /work/PDLINPUT_SOURCE -C /input/nanoaod.yaml --mc -O "/out/$1"
        /build/delphi-raw-nanoaod/delphi-raw-nanoaod -P /work/PDLINPUT_SOURCE -O "/out/$2"
        ' bash "$native" "nanoaod_raw_$source.root"
done

# Fortran can return zero on input failure: check actual physics events.
uv run python - "$final" "$events" <<'PY'
from pathlib import Path
import sys
import uproot

directory, expected = Path(sys.argv[1]), int(sys.argv[2])
for name in ("nanoaod.root", "nanoaod_raw_sdst.root", "nanoaod_fadana.root", "nanoaod_raw_fadana.root"):
    with uproot.open(directory / name) as root:
        tree = root["Events"]
        count = tree.num_entries
        # Raw run records can retain positive event numbers but have no MC event.
        if "_raw_" in name:
            count = int((tree["Event_isMC"].array(library="np") == 1).sum())
    if count != expected:
        raise RuntimeError(f"{name}: expected {expected}, found {count}")
    print(f"{name}: {count} physics events")
PY
touch "$work/.done"
