# DELPHI Analysis

Small analysis workspace for the merged SDST, RAWSDST, and RAWFADANA NanoAOD workflow.

## Setup

```bash
micromamba create -y -f environment.yml
source setup.sh
```

## Workflow

The scripts use `SAMPLE_SET = "100kTest"` by default. Edit the constants near the top of
each script when switching sample sets or output locations.

```bash
./scripts/filter-chunks.py
./scripts/hadd-chunks.py
./scripts/summary-branches.py
./scripts/plot-branches-all.py
./scripts/plot-branches-compare.py
./scripts/plot-gen-compare.py
./scripts/plot-gen-reco-track-match.py
./scripts/plot-gen-reco-track-diff.py
./scripts/plot-gen-check.py
```

The output layout is:

```text
output/chunks_raw/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
output/chunks/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
output/dataset/<sample-set>/<sample>/<nanoaod-file>.root
data/check/<sample-set>/
data/branch/<sample>.json
plots/<sample-set>/<plot-name>/
```

`filter-chunks.py` keeps events with `nGenPart > 0`, removes duplicate
`(run, event, nGenPart)` keys per job, and writes a per-job summary to
`data/check/<sample-set>/event-filter.csv`. `hadd-chunks.py` then combines the filtered
chunks into `output/dataset/<sample-set>/` without ROOT `hadd`.

The branch JSON files store per-sample Events branch stats by source file. `plot-branches-all.py`
and `plot-branches-compare.py` use those stats to skip branches with no filled entries.
Plotting scripts write diagnostic text or CSV files under `data/check/<sample-set>/<plot-name>/`.
