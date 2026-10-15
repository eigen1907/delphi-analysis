# PHOTOS-FSR production

Run from `delphi-analysis` on this Apple Silicon Mac, with Docker running:

```bash
mkdir -p data/20261014_jongwon
nohup caffeinate -i bash sample/run.sh > data/20261014_jongwon/run.log 2>&1 &
```

`nohup` keeps the process running after closing the terminal.
`caffeinate -i` prevents idle sleep while the command runs.
Follow progress with `tail -f data/20261014_jongwon/run.log`.
Default: `--sample ZKK --ncpu 8`. Select another sample or concurrency:

```bash
bash sample/run.sh --sample Zmumu --ncpu 4
```

`--ncpu` is the maximum number of concurrent jobs. Jobs run in batches; each
job runs Gen → simulation → conversion in order. Setup is performed once.
Run one launcher at a time for this campaign. An active old `00_run.sh` is
detected and blocks a new run; changing scripts does not change that process.

## Execution order

| Script | Task |
| --- | --- |
| `run.sh` | Run setup, then 20 jobs for the selected sample (8 at a time by default) |
| `01_setup.sh` | Build Docker/converter and prepare conditions |
| `02_generate.sh` | Generate PYTHIA + PHOTOS input, including the EOF buffer |
| `03_simulate.sh` | Run unchanged DELSIM → DELANA → SDST |
| `04_convert.sh` | Write native/raw ROOT and check event counts |

`config.env` contains host paths and image names. `Dockerfile` defines the
container's packages and environment. `01_setup.sh` builds the image, prepares
conditions and builds the converter.
Sources in `../delphisimulation` and its `delphi-nanoaod` checkout are read-only.
Builds and conditions stay in `.cache/sample/`.

## Event counts

Every sample uses **20 jobs × 5,000 simulated events = 100,000 events**.
There are no sample-specific event counts. A run processes one selected sample.
Gen writes `ceil(5000 / 0.9) = 5556` records/job as an EOF buffer; DELSIM
processes exactly 5,000. Use the detector file's Gen record for the same events.

PYTHIA seeds are `1907–2006`, DELSIM runs `3000–3099`, unique per job.
These remain tied to sample/job, independent of execution order or `--ncpu`.
Completed jobs get `production/job_N/.done` only after ROOT validation.
Rerunning `run.sh` with the same sample, event count and sources skips these
jobs. Incomplete jobs and old 4,500-event jobs are regenerated. If a job fails,
its active batch finishes before the launcher stops. Interrupting the launcher
also waits for its active workers. Logs and Gen/DST inputs are retained.

## Output

```text
data/20261014_jongwon/
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

After editing only the converter, rebuild once and reconvert the saved jobs:

```bash
bash sample/01_setup.sh > data/20261014_jongwon/setup.log 2>&1
for job in data/20261014_jongwon/*/production/job_*; do
    bash sample/04_convert.sh "$job" > "$job/conversion.log" 2>&1
done
```

## Physics notes

Four cards match the upstream PHOTOS settings. The original Florian tau
PHOTOS card was unavailable; the tau card uses the same radiation settings
and PYTHIA tau decays. This is a new comparable production.
Generator `Ecm=91.187 GeV` and DELSIM `EBEAM=45.625 GeV` preserve the existing
small energy difference. STIC truth links remain unresolved; raw full-DST
EM parsing remains incomplete (use raw SDST for that study).
