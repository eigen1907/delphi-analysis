# DELPHI Analysis

ROOT validation and event classification for DELPHI NanoAOD.

## Environment

Python 3.12 and the dependencies are defined in `pyproject.toml` and `uv.lock`:

```bash
uv sync --locked
```

Run commands from the repository root with `uv run python`. The reusable Python
packages are installed from `src/`; no `PYTHONPATH` changes are needed. ROOT is optional
for dataset preparation with `--backend hadd`; the default uses uproot.

## Layout

```text
scripts/
  02_sample_preparation/
    prepare_dataset.py          prepare chunks, merge and build feature datasets
  03_sample_validation/
    00_overall/
    01_track/
    02_gamma/
    03_rich/
  04_event_classification/
    01_model_train.py
    02_model_evaluate.py
    03_model_apply.py
src/
  common/                       paths, sample names and shared analysis helpers
  sample_preparation/           ROOT filtering and merging
  sample_validation/            overall, track, gamma and RICH implementations
  event_classification/         features, job splits, training, evaluation, inference
config/pid/                     feature and hyperparameter JSON configurations
docs/                           branch audits and converter notes
data/                           external and generated source samples
output/
  dataset/                      prepared ROOT files, feature datasets, predictions
  plots/                        validation and classifier figures
  models/                       models, model metadata and evaluation metrics
```

The numbered directories and files are execution entry points. Python packages
under `src/` use ordinary import names without numeric prefixes. Generated data,
outputs, container images and caches are excluded from Git.

## 1. Prepare datasets

One entry point provides the three preparation stages:

```bash
uv run python scripts/02_sample_preparation/prepare_dataset.py \
    -i data/20261015_jongwon -o output/dataset/20261015_jongwon
```

The default `--stage all` writes:

```text
output/dataset/20261015_jongwon/
  chunks/<sample>/final_root/job_N/   filtered per-job ROOT files
  merged/<sample>/                   merged analysis ROOT files
  classification/pid/                train.root, val.root, test.root, metadata.json
```

Use `--stage chunks`, `--stage merge` or `--stage classification` to run one
stage. The latter two accept the prepared chunks directory as input:

```bash
uv run python scripts/02_sample_preparation/prepare_dataset.py --stage merge \
    -i output/dataset/20261015_jongwon/chunks -o output/dataset/20261015_jongwon
uv run python scripts/02_sample_preparation/prepare_dataset.py --stage classification \
    -i output/dataset/20261015_jongwon/chunks -o output/dataset/20261015_jongwon \
    --feature-config config/pid/features/minimal.json
```

Filtering retains valid generator events and removes duplicated event keys.
The filtering summary is written under `output/dataset/checks/` unless `--check`
is specified. Merging is for analysis; classification reads the per-job chunks
and splits whole jobs into train/validation/test, preserving their independence.
The default split is 60/20/20 with seed 1907. `--samples` restricts filtering and
merging; classification requires all four classifier channels.

Feature configurations are `minimal`, `detector`, `pid` and `combined`. Their
names select the directory under `classification/`. Feature contents and order
are preserved in the dataset metadata and saved model.

## 2. Validate samples

Overall branch, generator and reconstructed distributions read merged data:

```bash
uv run python scripts/03_sample_validation/00_overall/00_branch_distribution.py \
    -i output/dataset/20261015_jongwon/merged -o output/plots/20261015_jongwon
```

Track, photon and RICH entry points accept source campaign directories with
per-job ROOT files. Each `plot.py` runs its stages while reading each sample once:

```bash
uv run python scripts/03_sample_validation/01_track/plot.py \
    -i data/20261015_jongwon -o output/plots
uv run python scripts/03_sample_validation/02_gamma/plot.py \
    -i data/20261015_jongwon -o output/plots
uv run python scripts/03_sample_validation/03_rich/plot.py \
    -i data/20261015_jongwon -o output/plots
```

The numbered entry points within each directory run a single stage. Their CLI
options and existing stage names/plot subdirectories are retained. Use `--help`
for sample selection and plotting options. Shared diagnostic tables are written
under `output/dataset/checks/`.

## 3. Train, evaluate and apply a classifier

The classifier uses four labels: `Zee`, `Zmumu`, `ZKK` and `Zpipi`. Tau samples
remain available for generation and validation. Train uses only train/validation;
validation is retained for early stopping. Evaluation loads the saved model and
reports train/validation/test metrics separately.

```bash
uv run python scripts/04_event_classification/01_model_train.py \
    -i output/dataset/20261015_jongwon/classification/pid \
    -o output/models/20261015_jongwon_pid_standard \
    --hyperparameters config/pid/hyperparameters/standard.json
uv run python scripts/04_event_classification/02_model_evaluate.py \
    -i output/dataset/20261015_jongwon/classification/pid \
    -m output/models/20261015_jongwon_pid_standard/bdt.joblib \
    -o output/plots/event_classification/20261015_jongwon_pid_standard
uv run python scripts/04_event_classification/03_model_apply.py \
    -i output/dataset/20261015_jongwon/merged \
    -m output/models/20261015_jongwon_pid_standard/bdt.joblib \
    -o output/dataset/20261015_jongwon/predictions/20261015_jongwon_pid_standard
```

Training writes `bdt.joblib`. Evaluation writes `metrics.json` and
`feature_importance.csv` beside the model and diagnostic PNGs to its plot output.
Apply writes prediction ROOT files without modifying inputs. Existing model
bundles remain readable. Match dataset features and feature configuration to the
model when evaluating or applying it.

## Data conventions

Source samples in `data/` are retained unchanged. Filtering and merging write
new files under `output/dataset/`; generation intermediate files stay with their
source jobs. A structure change does not recalibrate or repair physics content.

The existing matching algorithms, detector selections and numbered plot stages
are retained. DELPHI `Part_pdgId`/`SimPart_pdgId` use DELPHI mass codes, not PDG
codes. Raw run/setup rows are not physics events; production validation counts
`Event_isMC == 1`. RICH and truth-link limitations are recorded in the
[branch usage summary](docs/branch_usage_summary.md),
[raw SDST audit](docs/nanoaod_raw_sdst_branch_audit.md) and
[converter validation notes](docs/gamma_rich_converter_fixes.md).
