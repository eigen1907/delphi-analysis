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
        |          -> scripts/plot/*
        |               -> plots/<sample-set>
        |
        +-- scripts/bdt/prepare.py
              -> data/<study>/ml/<sample-set>_<feature-set>
                   -> scripts/bdt/train.py
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

Generated contents of `data/` and `plots/` are ignored by Git. Their `.gitkeep`
files preserve the directories. Study recipes in `runs/` are versioned.

Run the recipes from the repository root:

```bash
bash runs/prepare_data.sh
bash runs/plot_checks.sh
bash runs/plot_isr_photons.sh
bash runs/train_bdt_pid_standard.sh
bash runs/apply_bdt_pid_standard.sh
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

The ISR recipe reads `data/20260828_florian` and uses only its five
`photosFSR` samples: `Zee`, `Zmumu`, `Ztautau`, `ZKK`, and `Zpipi`.
The `stable_gen_*` plots include all gen particles with `GenPart_pdgId == 22`
and `GenPart_status == 1`, with no energy cut. Gen energy comes from
`GenPart_vector.fCoordinates.fT`. The `reco_*` plots use the
`Photon_fourMomentum.fCoordinates` branch: energy is `fT`, and the angles
come from `fX`, `fY`, and `fZ`. Both populations require a valid momentum
direction. The `isr_*` plots select stable gen photons with
$E_\gamma \geq 0.1$ GeV whose first non-photon ancestor is an electron
outside the hard-parent ancestry. The hard parent is the common parent of
the selected final-state pair. The `isr_matched_reco_*` plots show the reco
side of matches to those selected gen ISR photons. Each population has
multiplicity per event, $E_\gamma$, $\cos\theta_\gamma$, and $\phi_\gamma$
plots. The FSR implementation is parked in `python/plot_fsr_photons.py`; it has
no default recipe and does not produce plots during the ISR run.

Gen ISR photons are matched one to one to reco `Photon` candidates with an
opening angle below 0.05 rad. A populated `Photon`→`Part`→`SimPart`→`GenPart`
link that disagrees with an angular pair vetoes it; an absent link permits the
angular match. ISR and FSR candidates share the one-to-one matching step, so
one reco photon cannot be counted twice. Efficiency is matched gen ISR photons
divided by selected gen ISR photons, binned in gen energy or gen
$\cos\theta_\gamma$. No fiducial-angle cut is applied to the denominator.
Matching checks also show the nearest
opening angle before the matching cut and the reco-to-gen energy ratio for
matched photons.

The event energy plots compare the total selected gen ISR energy
$E_{\mathrm{ISR}}^{\mathrm{gen}}$ with the total energy of reco photons matched
to ISR, $E_{\mathrm{ISR}}^{\mathrm{reco}}$. In each bin of total gen ISR energy,
the energy recovery plot shows
$\sum E_{\mathrm{ISR}}^{\mathrm{reco}}/\sum E_{\mathrm{ISR}}^{\mathrm{gen}}$.
The radiated, recovered, and difference distributions show
$E_{\mathrm{ISR}}^{\mathrm{gen}}/\sqrt{s}$,
$E_{\mathrm{ISR}}^{\mathrm{reco}}/\sqrt{s}$, and
$(E_{\mathrm{ISR}}^{\mathrm{gen}}-E_{\mathrm{ISR}}^{\mathrm{reco}})/\sqrt{s}$
per event, using `Event_cmEnergy` as $\sqrt{s}$ (91.25 GeV in these samples).
The difference can be negative when reco energy exceeds gen energy; these
fractions do not directly determine the reconstructed collision energy. A
grouped bar plot shows the mean of each fraction per event in percent, making
the overall scale visible despite the large zero-event peak in the distributions.

Plots are grouped under `plots/20260828_florian/isr_photons/`:

```text
01_gen/         stable gen and selected ISR distributions
02_reco/        all reco Photon distributions
03_matching/    matched reco distributions and matching/energy checks
04_efficiency/  ISR count efficiency and energy recovery relative to gen ISR and sqrt(s)
```

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
