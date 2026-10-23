#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.env"

sample=ZKK
ncpu=8
events_per_job=5000
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) sample=${2:?--sample requires a sample name}; shift 2 ;;
        --ncpu) ncpu=${2:?--ncpu requires a positive integer}; shift 2 ;;
        -h|--help)
            echo "Usage: bash scripts/01_sample_generation/run.sh [--sample all|ZKK|Zpipi|Zmumu|Zee|Ztautau] [--ncpu 8]"
            echo "Defaults: ZKK, 8 concurrent jobs, 5000 events/job, 20 jobs."
            exit 0 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$sample" in
    all) samples=(ZKK Zpipi Zmumu Zee Ztautau) ;;
    ZKK|Zpipi|Zmumu|Zee|Ztautau) samples=("$sample") ;;
    *) echo "Unknown sample: $sample" >&2; exit 2 ;;
esac
if [[ ! "$ncpu" =~ ^[1-9][0-9]*$ ]]; then
    echo "--ncpu must be a positive integer" >&2
    exit 2
fi

# The old sequential runner does not use a lock and may still own these jobs.
if pgrep -f '[s]ample/00_run.sh' > /dev/null; then
    echo "00_run.sh is still running; finish or stop it before starting run.sh." >&2
    exit 1
fi
mkdir -p "$OUTPUT"
lock="$OUTPUT/.run.lock"
if ! mkdir "$lock"; then
    echo "Another run holds $lock; check for an active run before removing a stale lock." >&2
    exit 1
fi
# Wait for our workers before releasing the shared build/output directory.
trap 'wait; rmdir "$lock"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

echo "$(date): setup (see $OUTPUT/setup.log)"
if ! bash "$SAMPLE/01_setup.sh" > "$OUTPUT/setup.log" 2>&1; then
    echo "Setup failed; see $OUTPUT/setup.log" >&2
    exit 1
fi

run_job() {
    local sample=$1
    local offset=$2
    local job=$3
    local seed=$((1907 + offset + job))
    local nrun=$((3000 + offset + job))
    local work="$OUTPUT/20261014_100kTest_${sample}_photosFSR/production/job_$job"
    local log="$work/generator.log"
    trap 'echo "$(date): $sample job_$job failed; see $log" >&2' ERR

    mkdir -p "$work"
    printf '%s\n' "$events_per_job" > "$work/events.txt"
    echo "$(date): $sample job_$job generating ($events_per_job events; see $log)"
    bash "$SAMPLE/02_generate.sh" "$work" "$sample" "$seed" > "$log" 2>&1
    log="$work/simulation.log"
    echo "$(date): $sample job_$job simulating (see $log)"
    bash "$SAMPLE/03_simulate.sh" "$work" "$nrun" > "$log" 2>&1
    log="$work/conversion.log"
    echo "$(date): $sample job_$job converting (see $log)"
    bash "$SAMPLE/04_convert.sh" "$work" > "$log" 2>&1
    echo "$(date): $sample job_$job complete"
}

wait_batch() {
    local pid
    local failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then failed=1; fi
    done
    pids=()
    if [[ "$failed" == 1 ]]; then
        echo "A job failed; completed jobs are saved. Check the job log above." >&2
    fi
    return "$failed"
}

# All samples share the same limit on concurrent jobs.
pids=()
total_jobs=$((20 * ${#samples[@]}))
echo "$(date): $sample, $total_jobs jobs, up to $ncpu concurrent jobs"
for channel in "${samples[@]}"; do
    case "$channel" in
        ZKK) offset=0 ;;
        Zpipi) offset=20 ;;
        Zmumu) offset=40 ;;
        Zee) offset=60 ;;
        Ztautau) offset=80 ;;
    esac
    for job in {0..19}; do
        work="$OUTPUT/20261014_100kTest_${channel}_photosFSR/production/job_$job"
        if [[ -f "$work/.done" ]] && [[ "$(cat "$work/events.txt")" == "$events_per_job" ]]; then
            echo "$channel job_$job: already complete"
            continue
        fi
        run_job "$channel" "$offset" "$job" &
        pids+=("$!")
        if [[ ${#pids[@]} -ge "$ncpu" ]]; then wait_batch; fi
    done
done
if [[ ${#pids[@]} -gt 0 ]]; then wait_batch; fi
echo "$(date): $sample complete, $((total_jobs * events_per_job)) events: $OUTPUT"
