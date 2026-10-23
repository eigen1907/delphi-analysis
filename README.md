# DELPHI Analysis

ROOT validation and BDT event classification for DELPHI NanoAOD.

## Setup

Use Python 3.12 and run commands from the repository root:

```bash
uv sync --locked
```

```text
scripts/
  02_sample_preparation/      paired train/val/test datasets
  03_sample_validation/       overall, track, gamma and RICH plots
  04_event_classification/    BDT training, evaluation and inference
src/                         reusable Python code
config/pid/                  feature sets and hyperparameters
docs/                        data definitions and branch checks
data/                        source samples
output/{dataset,plots,models}/
```

Bash recipes contain defaults; append arguments to override them. Packages are
installed from `src/`, so no `PYTHONPATH` setting is needed. Data, outputs
and caches are excluded from Git.

## Prepare datasets

Pair native/raw events and write one train/val/test dataset per sample:

```bash
bash scripts/02_sample_preparation/prepare_dataset.sh \
    --input data/20261015_jongwon --output output/dataset/20261015_jongwon \
    --counts 60000 15000 15000 --seed 1907
```

`--counts` means TRAIN VAL TEST. Sources are preserved and existing datasets
are not overwritten. See [data rules](docs/data.md) for the Florian commands
and its 14,999-event Zmumu test split.

```text
output/dataset/<campaign>/
  metadata.json
  <sample>_photosFSR/{train,val,test}.root
    analysisLevel/Events
    raw_sdst/Events
```

The two trees contain the same events in the same order, keyed by
`(job_id, run, event)`. Preparation uses uproot and needs no ROOT installation.

## Validate samples

```bash
bash scripts/03_sample_validation/plot.sh --input output/dataset/20261015_jongwon
bash scripts/03_sample_validation/plot.sh --input output/dataset/20260828_florian
```

Plots use `train` by default; pass `--split val` or `--split test` to change it.
Each study's `plot.sh` runs all stages; numbered recipes run one stage.
Outputs follow the script names:

```text
output/plots/<campaign>/03_sample_validation/<split>/
  {00_overall,01_track,02_gamma,03_rich}/<numbered_stage>/
```

Both campaigns use a 1° gamma matching cone. See [photon definitions](docs/photons.md)
for selections and efficiency denominators. Old Florian RICH associations
support storage checks only; see [converter changes](docs/converter.md).

## Train and evaluate

Choose a campaign and one of `detector`, `detector_pid` (default),
`detector_rich` or `detector_pid_rich`. All use the `default` hyperparameters:

```bash
export CAMPAIGN=20261015_jongwon FEATURE_SET=detector_pid
bash scripts/04_event_classification/01_model_train.sh
bash scripts/04_event_classification/02_model_evaluate.sh
bash scripts/04_event_classification/03_model_apply.sh
```

Training fits `train` and uses `val` for early stopping. Use `test` for the final comparison.
See [BDT inputs](docs/bdt.md) for features, classes and missing-value rules.

| Output | Location |
| --- | --- |
| Models and metrics | `output/models/<campaign>_<feature>_default/` |
| Evaluation plots | `output/plots/<campaign>/04_event_classification/02_model_evaluate/<feature>_default/` |
| Predictions | `output/dataset/<campaign>/predictions/<feature>_default/` |

Experiment summaries are stored in `output/models/results.csv`.
See the [documentation index](docs/README.md) for branch checks and known limits.
