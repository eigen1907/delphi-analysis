#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
project="$PWD"
cache="$project/.cache/short94_c2"

# Reuse the prepared Florian-baseline converter, without the Gamma/RICH fixes.
convert() {
    local directory=$1
    local name=${directory##*/}
    local work="$cache/work/$name"
    local log="$directory/conversion.log"
    mkdir -p "$work"
    printf 'FILE=/dataset/%s.al\n' "$name" > "$work/PDLINPUT_SOURCE"
    echo "$(date): $name converting"
    {
        cat "$cache/provenance.txt"
        echo "Original upstream configuration: --data, select_hadrons=false, IFLRNQ=0"
        echo "# RAW_RECORD lines identify stored record types; no physics selection added."
        docker run --init --rm --platform linux/amd64 --name "short94-$name" \
            -v "$cache/build:/build:ro" -v "$cache/src:/src:ro" \
            -v "$project/.cache/sample/condition-data:/eos/opendata/delphi/condition-data:ro" \
            -v "$directory:/dataset" -v "$work:/work" delphi-analysis:local -c '
            set -e
            /build/delphi-nanoaod/delphi-nanoaod -P /work/PDLINPUT_SOURCE \
                -C /src/config/delphi-nanoaod.yaml --data -O /dataset/nanoaod.root
            /build/delphi-raw-nanoaod/delphi-raw-nanoaod -P /work/PDLINPUT_SOURCE \
                -O /dataset/nanoaod_raw_sdst.root
            '
    } > "$log" 2>&1
    uv run python "$cache/validate.py" "$directory" >> "$log" 2>&1
    tail -2 "$log"
    echo "$(date): $name complete"
}

pids=()
for input in "$project"/data/20260828_florian/20261026_short94_c2/*/*.al; do
    directory=${input%/*}
    if [[ -s "$directory/nanoaod.root" && -s "$directory/nanoaod_raw_sdst.root" ]] &&
       grep -q "Native/raw physics keys: exact 1:1" "$directory/conversion.log"; then
        echo "${directory##*/}: already validated"
        continue
    fi
    convert "$directory" &
    pids+=("$!")
    if [[ ${#pids[@]} -eq 3 ]]; then
        for pid in "${pids[@]}"; do wait "$pid"; done
        pids=()
    fi
done
if [[ ${#pids[@]} -gt 0 ]]; then
    for pid in "${pids[@]}"; do wait "$pid"; done
fi
