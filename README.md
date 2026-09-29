# DELPHI Analysis

Compact analysis workspace for DELPHI SDST, RAW-SDST, and RAW-FADANA
NanoAOD samples.

## Environment

The Python dependencies are defined in `pyproject.toml` and pinned in `uv.lock`.
Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run:

```bash
uv sync --locked
```

On macOS, install the OpenMP runtime with `brew install libomp` before using
XGBoost. ROOT is only needed for the optional `--backend hadd` merge mode;
the default uproot backend works with the uv environment alone.

## Workflow

```text
data/<study>/raw/<sample-set>
  |
  +-- scripts/data/prepare-chunks.py
  |     -> data/<study>/chunk/<sample-set>
  |
  +-- data/<study>/chunk/<sample-set>
        |
        +-- scripts/data/build-dataset.py
        |     -> data/<study>/dataset/<sample-set>
        |          -> scripts/<study>/*
        |               -> plots/<sample-set>
        |
        +-- scripts/pid/prepare.py
              -> data/<study>/ml/<sample-set>_<feature-set>
                   -> scripts/pid/train.py
                        -> plots/bdt/<sample-set>_<feature-set>_<profile>
```

The repository keeps code, recipes, and generated files in separate locations:

| Directory | Purpose |
| --- | --- |
| `python/` | Reusable analysis and feature construction code |
| `scripts/` | Command-line programs for data, plots, and BDTs |
| `config/` | Feature sets and hyperparameter profiles |
| `runs/` | Reproducible study commands |
| `data/` | Input data, intermediate files, and ML datasets |
| `plots/` | Generated figures and BDT run outputs |

Code and recipes are grouped by study with matching directory names:

```text
python/delphi_analysis/
  plot_utils.py       shared plotting helpers
  checks/             branch, generator, and reco validation
  tracking/           gen–reco track matching
  rich/               RICH inspection
  pid/                BDT preparation, training, and application
  photons/            photon reconstruction and ISR study
scripts/{checks,tracking,rich,pid,photons,data}/
runs/{checks,tracking,rich,pid,photons,data}/
config/pid/{features,hyperparameters}/
```

Generated contents of `data/` and `plots/` are ignored by Git. Their `.gitkeep`
files preserve the directories. Study recipes in `runs/` are versioned.

Run the recipes from the repository root:

```bash
bash runs/data/prepare_data.sh
bash runs/checks/plot_checks.sh
bash runs/tracking/plot_track_matching.sh
bash runs/rich/plot_rich.sh
bash runs/photons/plot_isr_photons.sh
bash runs/pid/train_bdt_pid_standard.sh
bash runs/pid/apply_bdt_pid_standard.sh
```

The data preparation, checks, and BDT recipes target `20260606_100kTest`
under `data/202606xx_jongwon/`. The photon recipe targets the Florian samples.
Edit the sample and data path in a recipe for another study.

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
isr-photons.py                  generated ISR photons and reconstruction efficiency
```

Plotting diagnostics and CSV summaries are written below
`data/check/<sample-set>/`.

## Photon and ISR study

```bash
bash runs/photons/plot_isr_photons.sh
```

The recipe reads the five Florian `photosFSR` samples: `Zee`, `Zmumu`,
`Ztautau`, `ZKK`, and `Zpipi`. The code lives in
`python/delphi_analysis/photons/`: `data.py` reads and classifies photons,
`plot_isr_photons.py` makes 1D figures, and `plot_2d.py` makes energy–angle
maps. The FSR-specific script is reserved for a later study.

All stable gen photons and all reco photons enter the study. There is no gen
energy threshold or fiducial-angle cut. ISR is divided into beam-collinear
photons and the remaining non-collinear photons. Stable gen photons are
partitioned into non-collinear ISR, beam-collinear ISR, FSR from the hard pair,
and all other photons. The exact branch and ancestry rules are in `data.py`.
Four-vectors must have finite, positive energy and nonzero momentum; invalid
records stop the run. None were found in these five samples.

Energy, cos(theta), and phi use component stacks with the same event
normalization. An exactly beam-directed photon has undefined phi; the stack
**displays** it in the phi=0 bin, as the figure title states. That location is
not a physical angle measurement. Multiplicity plots overlay the total and
component distributions because their histogram heights cannot be stacked to
obtain the total multiplicity distribution.

All stable gen and reco photons compete in one closest-angle-first one-to-one
match. The opening angle must be below 0.05 rad. A known truth association to
another gen particle vetoes the pair; an absent association permits angular
matching. No energy compatibility cut is applied. The `truth_linked/` companion
figures contain only matched pairs whose stored association resolves to the
same gen photon. Their denominators remain inclusive, so the companion figures
also depend on truth-association completeness.

Inclusive photon and ISR efficiencies use matched gen photons over all gen
photons of that population, with gen coordinates in both numerator and
denominator. The 1D plots show 68.27% Clopper–Pearson intervals. The 2D maps
include the denominator and interval width, with empty denominator bins shown
gray. These efficiencies include detector acceptance. The beam-collinear
records are included in the denominators.

The inclusive angular match can associate very soft beam-collinear records
with unrelated forward reco photons. In these samples, 225 such pairs have no
truth association. Therefore the angular-associated reco energy is a matching
diagnostic, not a validated physical ISR energy recovery. The `truth_linked/`
companion results expose this difference. Total gen ISR energy always includes
all ISR photons, including the beam-collinear component. Energy fractions use
the per-event `Event_cmEnergy`; reco energies can exceed gen energies and are
not clipped.

Plots are grouped under `plots/20260828_florian/isr_photons/`:

```text
01_gen/         stable gen and ISR distributions, component stacks, gen 2D maps
02_reco/        reco photon distributions and 2D map
03_matching/    matched gen/reco distributions and matching 2D maps
04_efficiency/  1D/2D photon and ISR efficiencies, energy accounting,
                and truth_linked/ companion results
```

## BDT Classification

The classifier uses `xgboost.XGBClassifier`. Four classes are assigned fixed labels:
`Zee`, `Zmumu`, `ZKK`, and `Zpipi`. Jobs are shuffled with a fixed seed and split
independently per class; the default train/validation/test fractions are 60/20/20.

Feature configurations are stored separately:

```text
config/pid/features/minimal.json   focused tracking, calorimeter, and vertex detector baseline
config/pid/features/detector.json  all combined inputs except direct DELPHI PID decisions
config/pid/features/pid.json       tracking, vertex, and DELPHI PID outputs
config/pid/features/combined.json  detector inputs plus MuidRaw, ElidRaw, and HaidRaw PID outputs
```

The JSON keys retain the ROOT collection names. Vector branches are ordered,
truncated to ten objects, and padded with `NaN`. Track-associated collections are
matched to the leading `TracRaw` entries through their association indices. The raw
15-element track covariance is retained for the four leading tracks. No sums, means,
ratios, pair variables, or reconstructed momentum features are calculated.

Hyperparameter profiles are stored in:

```text
config/pid/hyperparameters/light.json
config/pid/hyperparameters/standard.json
config/pid/hyperparameters/heavy.json
```

`prepare.py` writes `train.root`, `val.root`, `test.root`, and `metadata.json`.
The metadata records the class mapping, selected jobs, event counts, split fractions,
feature configuration, and exact expanded feature list.

`train.py` writes the fitted model, metrics, normalized validation and test
confusion matrices, and XGBoost gain feature importance under `--output`.
`apply.py` writes event identifiers, the predicted class, and per-class
probabilities without modifying the input ROOT file.

Generator truth, event identifiers, and MC-only metadata are not model inputs.

## Source Integrity

Files below `data/<study>/raw/` are treated as the source of truth. Filtering,
merging, plotting, and BDT preparation do not recalibrate or repair branch
contents. When the NanoAOD schema or converter changes, regenerate the raw
data before rebuilding the downstream products.

In particular, RICH measurements stored in `HaidRaw_*` should only be used after the
upstream `QGRIC/KGRIC` and `QLRIC/KLRIC` mappings and output types have been
validated. The current PID and combined feature configs therefore keep the DELPHI
PID decisions but exclude the unverified RICH ring summaries and quality word.
