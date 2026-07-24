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
./scripts/filter-chunks.py -i data/raw/<sample-set> -o data/chunk/<sample-set>
./scripts/hadd-chunks.py -i data/chunk/<sample-set> -o data/dataset/<sample-set>
```

Plotting scripts require `-i/--input`. If `-o/--output` is omitted, plots are written
under `plots/<input-directory-name>/`; for example, `-i /path/to/OpenData` writes to
`plots/OpenData/`.

```bash
./scripts/plot-branches-all.py -i data/dataset/<sample-set>
./scripts/plot-branches-compare.py -i data/dataset/<sample-set>
./scripts/plot-gen-compare.py -i data/dataset/<sample-set>
./scripts/plot-gen-reco-track-match-cut.py -i data/dataset/<sample-set>
./scripts/plot-gen-reco-track-match-result.py -i data/dataset/<sample-set>
./scripts/plot-gen-check.py -i data/dataset/<sample-set>
./scripts/plot-reco-check.py -i data/dataset/<sample-set>
```

Use `--samples sample_a sample_b` to process an explicit subset. Use `--data-root` or `--check`
when a specific script needs a custom diagnostic path.

The output layout is:

```text
data/raw/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
data/chunk/<sample-set>/<generated-sample>/final_root/job_<n>/<nanoaod-file>.root
data/dataset/<sample-set>/<sample>/<nanoaod-file>.root
data/ml/<sample-set>/{train,val,test}.root
data/check/<sample-set>/
plots/<sample-set>/<plot-name>/
```

`filter-chunks.py` keeps events with `nGenPart > 0`, removes duplicate
`(run, event, nGenPart)` keys per job, and writes a per-job summary to
`data/check/<sample-set>/event-filter.csv`. `hadd-chunks.py` then combines the filtered
chunks into `data/dataset/<sample-set>/` with ROOT `hadd` by default.

`plot-branches-all.py` and `plot-branches-compare.py` inspect the ROOT files directly
before plotting, so no separate branch summary JSON step is needed. Plotting scripts
write diagnostic text or CSV files under `data/check/<sample-set>/<plot-name>/`.

## BDT classification

The BDT input is built from the filtered `nanoaod_raw_sdst.root` chunks. Jobs are
shuffled with a fixed seed and split independently for each class, so events from one
job cannot appear in multiple splits.

```bash
./scripts/prepare-bdt-data.py \
  -i data/chunk/20260606_100kTest \
  -o data/ml/20260606_100kTest

./scripts/train-bdt.py \
  -i data/ml/20260606_100kTest \
  -o data/bdt/20260606_100kTest
```

The default split is 60% training, 20% validation, and 20% test by job for each
of `Zee`, `Zmumu`, `ZKK`, and `Zpipi`. The class labels, selected jobs, event
counts, and feature list are recorded in `metadata.json`.

BDT features are defined and extracted in `python/bdt_features.py`. The current
feature set uses the four highest-momentum tracks with matched PID and calorimeter
information, the four highest-energy EM and hadronic showers, and a small set of
event-level energy and momentum features. Generator truth and MC-only metadata are
excluded; event identifiers are stored only for tracing events and are not model
inputs.

Training writes the fitted model, metrics, normalized validation and test confusion
matrices, and validation permutation importance under the requested output directory.

Apply the trained model to another `nanoaod_raw_sdst.root` file with:

```bash
./scripts/apply-bdt.py \
  -i data/dataset/20260606_100kTest/Zee/nanoaod_raw_sdst.root \
  -m data/bdt/20260606_100kTest/bdt.joblib \
  -o data/bdt/20260606_100kTest/prediction_Zee.root
```

The prediction tree contains the run and event identifiers, `predicted_label`, and
one probability branch per class. The input ROOT file is not modified.
