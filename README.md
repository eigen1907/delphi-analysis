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

Stable gen photons have `GenPart_pdgId == 22` and `GenPart_status == 1`;
reco photons use `Photon_fourMomentum`. Gen ancestry divides stable photons
into Non-beam ISR, Beam ISR, FSR, and Others. Beam ISR has zero transverse
momentum and a direct parent from the incoming parentless electron or
positron. Non-beam ISR is the remaining ISR. The exact ancestry rules are in
`data.py`. There is no gen energy threshold, fiducial-angle cut, or matching
energy cut. Invalid photon four-vectors stop the run rather than being removed.

Stage 01 stacks the four stable-gen origin categories in the 1D energy,
cos(theta), and phi distributions. Stage 02 removes Beam ISR from the gen
stack and gen total. Both stages overlay the same set of **all** reco photon
candidates as black points with Poisson error bars. Gen and reco use the same
bins within each stage and are divided by the same number of events. This
compares populations; reco candidates need not match the gen photons
underneath them. Multiplicity shows the separate gen-origin event
distributions, the gen total, and reco points. The origin distributions cannot
be stacked into the total multiplicity distribution.
Photon-coordinate histograms show `Nγ / bin per event` (bin count divided by
the number of events), while multiplicity uses `Event fraction`. A photon with
undefined phi is displayed in the phi=0 bin solely for plotting.

Each stage also shows the event-by-event energy ratio
`R_E = (sum of all reco photon energies) / (sum of selected stable-gen photon energies)`.
Stage 01 includes Beam ISR in the gen sum; stage 02 excludes it. Positive
ratios use common logarithmic bins across all samples and both stages, and
defined events with no reco photons appear as a separate marker at `R_E = 0`.
The ratio axis is linear near zero and logarithmic above one. Events
with zero gen energy have an undefined ratio; their count is printed on each
plot and they are excluded from the ratio histogram. All fractions use the
full sample event count. The ratio can exceed one and is not a matched-photon
energy recovery or reconstruction efficiency.

The three matching methods are:

1. **Truth only:** Follow `Photon_partIdx → Part_simIdx → SimPart_genIdx` to a
   stable gen photon. No opening-angle requirement is applied.
2. **Angle only:** Match all stable gen photons against all reco photons by
   closest-first, one-to-one 3D opening angle below 0.03 rad, ignoring truth
   links.
3. **Truth + angle:** Assign direct truth pairs first. For unmatched stable
   gen photons, use the same angular rule with reco photons whose truth link is
   absent (`< 0`). A reco photon linked to a different gen particle is not an
   angular fallback candidate.

Efficiency plots use two denominators: all gen ISR photons (including Beam
ISR) and Non-beam ISR photons. In each case the numerator is the matched gen
ISR subset from each method, binned by gen energy or gen cos(theta). All
efficiencies include detector acceptance and have 68.27% Clopper–Pearson
intervals. The pooled opening-angle scan in `03_matching/` shows how the
direct-link and angle-only pair counts change with the cut. At 0.03 rad,
6,525/6,757 (96.6%) directly linked ISR pairs lie within the cut, and angle
matching associates 3 Beam ISR photons to reco candidates with no gen link.
At 0.05 rad those numbers are 6,737/6,757 (99.7%) and 225. The Beam no-link
scan curve is part of the total no-link curve; an absent link does not prove
an angular pair is incorrect. Efficiency plots use symmetric-log axes near
zero to show small values without excluding any photons.

Plots are grouped under `plots/20260828_florian/isr_photons/`. Each sample has
one axis per 1D figure. The gen and reco 2D energy–cos(theta) maps in stages
01 and 02 share a two-panel figure, common bins, and one color scale across
the five samples within each stage. Matched-gen 2D maps in stage 03 also use
common bins and per-method color scales across samples.

```text
01_gen_reco/<sample>/              gen/reco overlays, 2D map, event energy ratio
02_gen_reco_wo_beamISR/<sample>/  same with Beam ISR removed from gen
03_matching/<sample>/              matched gen ISR distributions and three 2D maps
03_matching/                       pooled ISR opening-angle cut scan
04_efficiency/<sample>/            all-ISR and Non-beam-ISR efficiencies vs energy and cos(theta)
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
