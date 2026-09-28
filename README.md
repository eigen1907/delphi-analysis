# DELPHI Analysis

Compact analysis workspace for DELPHI SDST, RAW-SDST, and RAW-FADANA
NanoAOD samples.

## Environment

The Python dependencies are defined in `pyproject.toml` and pinned in `uv.lock`.
Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run:

```bash
uv sync --locked
uv run --locked python scripts/plot/rich.py --help
```

On macOS, install the OpenMP runtime with `brew install libomp` before using
XGBoost. ROOT is only needed for the optional `--backend hadd` merge mode;
the default uproot backend works with the uv environment alone.

For an interactive shell, `source setup.sh` syncs the locked environment and
activates `.venv/`. Commands run with `uv run` do not need this step.

## Workflow

```text
data/raw/<sample-set>
  |
  +-- scripts/data/prepare-chunks.py
  |     -> data/chunk/<sample-set>
  |
  +-- data/chunk/<sample-set>
        |
        +-- scripts/data/build-dataset.py
        |     -> data/dataset/<sample-set>
        |          -> scripts/plot/*
        |               -> plots/<sample-set>
        |
        +-- scripts/bdt/prepare.py
              -> data/ml/<sample-set>_<feature-set>
                   -> scripts/bdt/train.py
                        -> plots/bdt/<sample-set>_<feature-set>_<profile>
```

Scripts are grouped by responsibility:

```text
scripts/data/  raw input filtering and dataset merging
scripts/plot/  dataset validation and diagnostic plots
scripts/bdt/   BDT preparation, training, and application
python/        shared implementation
config/bdt/    feature and hyperparameter configurations
```

## Data Preparation

`prepare-chunks.py` preserves the per-job directory structure. For trees with
generator information, it keeps events with `nGenPart > 0` and rejects duplicated
`(run, event, nGenPart)` keys. The filtering summary is written under
`data/check/<sample-set>/`.

`build-dataset.py` merges the prepared jobs independently for each sample and ROOT
file type. The uproot backend is used by default because ROOT's experimental RNTuple
merger may abort on these files. ROOT `hadd` remains available for compatible inputs.

The merged datasets are intended for validation plots and event-level inspection.
BDT preparation reads the prepared chunks directly so that complete jobs, rather
than individual events, can be assigned to train, validation, and test splits.

## Plotting

The plotting scripts inspect ROOT branches directly and write figures below
`plots/<sample-set>/`:

```text
branches-all.py                 every numeric branch by sample
branches-compare.py             common branches across NanoAOD sources
gen-check.py                    generator-level validation
gen-compare.py                  generator-level source comparison
reco-check.py                   reconstruction-level validation
rich.py                         RICH storage, dtype, and consistency study
gen-reco-track-match-cut.py     matching efficiency and cut scan
gen-reco-track-match-result.py  matched-track residuals
```

Plotting diagnostics and CSV summaries are written below
`data/check/<sample-set>/`.

## BDT Classification

The classifier uses `xgboost.XGBClassifier`. Four classes are assigned fixed labels:
`Zee`, `Zmumu`, `ZKK`, and `Zpipi`. Jobs are shuffled with a fixed seed and split
independently per class; the default train/validation/test fractions are 60/20/20.

Feature configurations are stored separately:

```text
config/bdt/features/minimal.json   focused tracking, calorimeter, and vertex detector baseline
config/bdt/features/detector.json  all combined inputs except direct DELPHI PID decisions
config/bdt/features/pid.json       tracking, vertex, and DELPHI PID outputs
config/bdt/features/combined.json  detector inputs plus MuidRaw, ElidRaw, and HaidRaw PID outputs
```

The JSON keys retain the ROOT collection names. Vector branches are ordered,
truncated to ten objects, and padded with `NaN`. Track-associated collections are
matched to the leading `TracRaw` entries through their association indices. The raw
15-element track covariance is retained for the four leading tracks. No sums, means,
ratios, pair variables, or reconstructed momentum features are calculated.

Hyperparameter profiles are stored in:

```text
config/bdt/hyperparameters/light.json
config/bdt/hyperparameters/standard.json
config/bdt/hyperparameters/heavy.json
```

`prepare.py` writes `train.root`, `val.root`, `test.root`, and `metadata.json`.
The metadata records the class mapping, selected jobs, event counts, split fractions,
feature configuration, and exact expanded feature list.

`train.py` writes the fitted model, metrics, normalized validation and test confusion
matrices, and XGBoost gain feature importance. `apply.py` writes event identifiers,
the predicted class, and per-class probabilities without modifying the input ROOT
file.

Generator truth, event identifiers, and MC-only metadata are not model inputs.

## Source Integrity

Files below `data/raw/` are treated as the source of truth. Filtering, merging,
plotting, and BDT preparation do not recalibrate or repair branch contents. When the
NanoAOD schema or converter changes, regenerate `data/raw/` before rebuilding the
downstream products.

In particular, RICH measurements stored in `HaidRaw_*` should only be used after the
upstream `QGRIC/KGRIC` and `QLRIC/KLRIC` mappings and output types have been
validated. The current PID and combined feature configs therefore keep the DELPHI
PID decisions but exclude the unverified RICH ring summaries and quality word.
