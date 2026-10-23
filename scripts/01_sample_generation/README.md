# Sample generation

Run from `delphi-analysis` on x86_64 Linux. One Apptainer image runs PYTHIA +
PHOTOS, DELSIM, DELANA and conversion. Set paths in [config.env](config.env).
Builds and conditions use `.cache/sample/`; simulation sources are read-only.

## Prepare the image

Build locally with Docker, using `linux/amd64` even on Apple Silicon:

```bash
docker buildx build --platform linux/amd64 --load \
    --build-context generator=../delphi-simulation/pythia8 \
    -t delphi-analysis:local scripts/01_sample_generation
docker save -o delphi-analysis.tar delphi-analysis:local
```

Adjust `--build-context` for a different checkout name. Save an existing amd64
image directly if suitable. Transfer the archive to the server, then run:

```bash
mkdir -p .cache/sample
apptainer build .cache/sample/delphi-analysis.sif docker-archive:delphi-analysis.tar
apptainer exec --cleanenv .cache/sample/delphi-analysis.sif /bin/bash -c \
    'set -e; command -v pythia8_generate; command -v runsim; command -v cmake'
git -C ../delphi-simulation submodule update --init --recursive
uv sync --locked
```

The server needs no Docker, fakeroot or CVMFS. Use `docker save` to preserve image
environment variables. First setup downloads about 3.2 GB of conditions.

## Run

All five samples, with up to 100 concurrent jobs:

```bash
mkdir -p data/20261015_jongwon
nohup bash scripts/01_sample_generation/run.sh --sample all --ncpu 100 > data/20261015_jongwon/run.log 2>&1 &
```

Follow `tail -f data/20261015_jongwon/run.log`; `nohup` survives logout.
Default: `--sample ZKK --ncpu 8`. Select another sample with `--sample Zmumu`.

`all` selects `ZKK`, `Zpipi`, `Zmumu`, `Zee`, `Ztautau`: **20 jobs × 5,000 events
per sample**, or 100 jobs / 500,000 events total. `--ncpu` limits concurrent jobs
across all samples, not CPU affinity. Jobs run in batches after one setup.
Run only one launcher per campaign; an active legacy `00_run.sh` blocks startup.

`run.sh` calls `01_setup.sh` once, then `02_generate.sh` → `03_simulate.sh` →
`04_convert.sh` per job. Setup builds the converter and prepares conditions.

Before full production, test ten events:

```bash
bash scripts/01_sample_generation/01_setup.sh
mkdir -p .cache/sample/smoke/production/job_0
printf '10\n' > .cache/sample/smoke/production/job_0/events.txt
bash scripts/01_sample_generation/02_generate.sh .cache/sample/smoke/production/job_0 ZKK 1907
bash scripts/01_sample_generation/03_simulate.sh .cache/sample/smoke/production/job_0 3000
bash scripts/01_sample_generation/04_convert.sh .cache/sample/smoke/production/job_0
```

## Output

Campaign files go to `data/20261015_jongwon/`: launcher/setup logs,
`nanoaod.patch` and `20261014_100kTest_<sample>_photosFSR/` for each sample.
Each sample has `production/job_0` … `job_19` and matching `final_root/job_N/`.
ROOT files: `nanoaod.root`, `nanoaod_raw_sdst.root`, `nanoaod_fadana.root`,
`nanoaod_raw_fadana.root`, `pythiastandalone.root`.

Each `production/job_N/` retains these intermediate files:

| File | Contents |
| --- | --- |
| `my_events.fadgen` | Generator input for DELSIM, including the EOF buffer |
| `simana.fadsim` | DELSIM output |
| `simana.fadana` | DELANA output |
| `simana.sdst` | Short DST |
| `generator.log`, `simulation.log`, `conversion.log` | Stage logs |

Raw ROOT includes run/setup rows; physics-event validation uses `Event_isMC == 1`.
The launcher log reports progress and detailed failure-log paths.

## Restart or reconvert

A job gets `.done` only after all four detector ROOT files pass event-count
validation. Reruns skip `.done` jobs with the same event count; incomplete jobs
and old 4,500-event jobs are regenerated, replacing their intermediate files.
Source changes are not detected: remove affected `.done` files to regenerate.
A failure or interrupt waits for active workers before the launcher exits.

After converter-only changes, rebuild and reconvert saved jobs:

```bash
bash scripts/01_sample_generation/01_setup.sh > data/20261015_jongwon/setup.log 2>&1
for job in data/20261015_jongwon/*/production/job_*; do
    bash scripts/01_sample_generation/04_convert.sh "$job" > "$job/conversion.log" 2>&1
done
```

Generator or Dockerfile changes require a local image rebuild and transfer;
remove the old SIF before converting the new archive. Converter-only changes
require setup, not an image rebuild.

## Production details

- Generation writes 5,556 records/job as an EOF buffer; DELSIM processes 5,000.
  Use the detector file's Gen record for paired events.
- PYTHIA seeds `1907–2006` and DELSIM runs `3000–3099` are unique, fixed by sample/job.
- Four cards match upstream PHOTOS settings. The original Florian tau card was
  unavailable; the tau card uses matching radiation settings and PYTHIA decays.
- Generator `Ecm=91.187 GeV` and DELSIM `EBEAM=45.625 GeV` retain their small
  energy difference. STIC truth links remain unresolved. Raw full-DST EM parsing
  is incomplete; use raw SDST for that study.
