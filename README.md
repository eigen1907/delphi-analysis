# DELPHI Analysis

Small analysis workspace for the merged SDST, RAWSDST, and RAWFADANA NanoAOD workflow.

## Setup

```bash
micromamba create -y -f environment.yml
source setup.sh
```

## Workflow

Chunk filtering and merging require explicit input and output roots. Sample names are
discovered from the input root, so different datasets can contain different samples.

```bash
./scripts/filter-chunks.py -i output/chunks_raw/<sample-set> -o output/chunks/<sample-set>
./scripts/hadd-chunks.py -i output/chunks/<sample-set> -o output/dataset/<sample-set>
```

Plotting scripts require `-i/--input`. If `-o/--output` is omitted, plots are written
under `plots/<input-directory-name>/`; for example, `-i /path/to/OpenData` writes to
`plots/OpenData/`.

```bash
./scripts/plot-branches-all.py -i output/dataset/<sample-set>
./scripts/plot-branches-compare.py -i output/dataset/<sample-set>
./scripts/plot-gen-compare.py -i output/dataset/<sample-set>
./scripts/plot-gen-reco-track-match-cut.py -i output/dataset/<sample-set>
./scripts/plot-gen-reco-track-match-result.py -i output/dataset/<sample-set>
./scripts/plot-gen-check.py -i output/dataset/<sample-set>
./scripts/plot-reco-check.py -i output/dataset/<sample-set>
```

Use `--samples sample_a sample_b` to process an explicit subset. Use `--data-root` or `--check`
when a specific script needs a custom diagnostic path.

The output layout is:

```text
output/chunks_raw/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
output/chunks/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
output/dataset/<sample-set>/<sample>/<nanoaod-file>.root
data/check/<sample-set>/
plots/<sample-set>/<plot-name>/
```

`filter-chunks.py` keeps events with `nGenPart > 0`, removes duplicate
`(run, event, nGenPart)` keys per job, and writes a per-job summary to
`data/check/<sample-set>/event-filter.csv`. `hadd-chunks.py` then combines the filtered
chunks into `output/dataset/<sample-set>/` without ROOT `hadd`.

`plot-branches-all.py` and `plot-branches-compare.py` inspect the ROOT files directly
before plotting, so no separate branch summary JSON step is needed. Plotting scripts
write diagnostic text or CSV files under `data/check/<sample-set>/<plot-name>/`.
