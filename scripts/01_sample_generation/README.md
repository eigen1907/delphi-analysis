# PHOTOS-FSR production

Run from `delphi-analysis` on an x86_64 Linux server. The generator, simulation
and converter use the same Apptainer SIF image.

Build the Docker image locally, then transfer it to this server and convert it
to SIF. Commands below run from `delphi-analysis` on each machine.

On the local machine, with Docker running:

```bash
docker buildx build --platform linux/amd64 --load \
    --build-context generator=../delphi-simulation/pythia8 \
    -t delphi-analysis:local scripts/01_sample_generation
docker save -o delphi-analysis.tar delphi-analysis:local
```

If the local simulation checkout has another name, adjust `--build-context`.
Use `linux/amd64` even on an Apple Silicon Mac, because the server is x86_64.
An existing `delphi-analysis:local` built from this Dockerfile for `linux/amd64`
can be saved directly without rebuilding.

Copy `delphi-analysis.tar` to the server's `delphi-analysis` directory, then run:

```bash
mkdir -p .cache/sample
apptainer build .cache/sample/delphi-analysis.sif docker-archive:delphi-analysis.tar
apptainer exec --cleanenv .cache/sample/delphi-analysis.sif /bin/bash -c \
    'set -e; command -v pythia8_generate; command -v runsim; command -v cmake'
git -C ../delphi-simulation submodule update --init --recursive
uv sync --locked
```

This conversion uses neither Docker nor fakeroot on the server. Package
installation and generator compilation have already finished in local Docker.
The server's `setgroups: Permission denied` failure occurs when running the
package-installation stage of a definition file, which this workflow avoids.
Use `docker save` to preserve the image's environment variables.

Run all five samples with up to 100 concurrent jobs:

```bash
mkdir -p data/20261015_jongwon
nohup bash scripts/01_sample_generation/run.sh --sample all --ncpu 100 > data/20261015_jongwon/run.log 2>&1 &
```

`nohup` keeps the process running after closing the terminal.
Follow progress with `tail -f data/20261015_jongwon/run.log`.
Default: `--sample ZKK --ncpu 8`. Select another sample or concurrency:

```bash
bash scripts/01_sample_generation/run.sh --sample Zmumu --ncpu 4
```

`--sample all` selects `ZKK`, `Zpipi`, `Zmumu`, `Zee` and `Ztautau`: 100 jobs
and 500,000 events in total. `--ncpu` is the maximum number of concurrent jobs
across all selected samples, rather than a CPU affinity setting. With
`--sample all --ncpu 100`, all 100 unfinished jobs can start together.
Jobs run in batches; each job runs Gen → simulation → conversion in order.
Setup is performed once.
Run one launcher at a time for this campaign. An active old `00_run.sh` is
detected and blocks a new run; changing scripts does not change that process.

## Execution order

| Script | Task |
| --- | --- |
| `run.sh` | Run setup, then 20 jobs per selected sample (8 at a time by default) |
| `01_setup.sh` | Check the SIF, build the converter and prepare conditions |
| `02_generate.sh` | Generate PYTHIA + PHOTOS input, including the EOF buffer |
| `03_simulate.sh` | Run unchanged DELSIM → DELANA → SDST |
| `04_convert.sh` | Write native/raw ROOT and check event counts |

`config.env` contains host paths and the SIF path. `Dockerfile` defines the
container's packages and environment. `01_setup.sh` checks that the SIF exists,
prepares conditions and builds the converter. After changing the generator or
`Dockerfile`, rebuild locally, transfer the new archive and convert it again
(remove the old SIF first). Converter source changes only require rerunning
setup.
Sources in `../delphi-simulation` and its `delphi-nanoaod` checkout are read-only.
Builds and conditions stay in `.cache/sample/`.
The first setup downloads about 3.2 GB of condition data; a CVMFS mount is not
required.

Before a full run, check all three stages with ten events:

```bash
bash scripts/01_sample_generation/01_setup.sh
mkdir -p .cache/sample/smoke/production/job_0
printf '10\n' > .cache/sample/smoke/production/job_0/events.txt
bash scripts/01_sample_generation/02_generate.sh .cache/sample/smoke/production/job_0 ZKK 1907
bash scripts/01_sample_generation/03_simulate.sh .cache/sample/smoke/production/job_0 3000
bash scripts/01_sample_generation/04_convert.sh .cache/sample/smoke/production/job_0
```

## Event counts

Every sample uses **20 jobs × 5,000 simulated events = 100,000 events**.
There are no sample-specific event counts. A run processes one sample or all five.
Gen writes `ceil(5000 / 0.9) = 5556` records/job as an EOF buffer; DELSIM
processes exactly 5,000. Use the detector file's Gen record for the same events.

PYTHIA seeds are `1907–2006`, DELSIM runs `3000–3099`, unique per job.
These remain tied to sample/job, independent of execution order or `--ncpu`.
Completed jobs get `production/job_N/.done` only after ROOT validation.
Rerunning `run.sh` with the same sample and event count skips these jobs.
Source changes are not checked; remove the affected jobs' `.done` files to
regenerate them after changing the generator, cards or converter.
Incomplete jobs and old 4,500-event jobs are regenerated. If a job fails,
its active batch finishes before the launcher stops. Interrupting the launcher
also waits for its active workers. Logs and Gen/DST inputs are retained.

## Output

```text
data/20261015_jongwon/
  run.log, setup.log, nanoaod.patch
  20261014_100kTest_ZKK_photosFSR/
    production/job_0/     # Gen, FADSIM, FADANA, SDST, configs and logs
    final_root/job_0/     # nanoaod.root, nanoaod_raw_sdst.root,
                         # nanoaod_fadana.root, nanoaod_raw_fadana.root,
                         # pythiastandalone.root
    ...                  # job_0 ... job_19
  ...                    # five samples
```

Raw ROOT also saves run/setup rows, sometimes retaining a previous positive
event number. Validation counts `Event_isMC == 1` as MC physics events.
The existing analysis reads this job-directory layout. `run.log` reports each
stage and, on failure, the path to its detailed log.

Intermediate files remain in each `production/job_N/` after ROOT conversion:

| File | Contents |
| --- | --- |
| `my_events.fadgen` | Generator input for DELSIM, including the EOF buffer |
| `simana.fadsim` | DELSIM output |
| `simana.fadana` | DELANA output |
| `simana.sdst` | Short DST |
| `generator.log`, `simulation.log`, `conversion.log` | Stage logs |

Restarting an incomplete job replaces its intermediate files. Completed jobs
with matching event counts are skipped.

After editing only the converter, rebuild once and reconvert the saved jobs:

```bash
bash scripts/01_sample_generation/01_setup.sh > data/20261015_jongwon/setup.log 2>&1
for job in data/20261015_jongwon/*/production/job_*; do
    bash scripts/01_sample_generation/04_convert.sh "$job" > "$job/conversion.log" 2>&1
done
```

## Physics notes

Four cards match the upstream PHOTOS settings. The original Florian tau
PHOTOS card was unavailable; the tau card uses the same radiation settings
and PYTHIA tau decays. This is a new comparable production.
Generator `Ecm=91.187 GeV` and DELSIM `EBEAM=45.625 GeV` preserve the existing
small energy difference. STIC truth links remain unresolved; raw full-DST
EM parsing remains incomplete (use raw SDST for that study).
