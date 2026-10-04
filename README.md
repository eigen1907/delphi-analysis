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
isr-photons.py                  photon topology and Gen–Sim–Reco truth efficiencies
```

Plotting diagnostics and CSV summaries are written below
`data/check/<sample-set>/`.

## Photon truth study

```bash
bash runs/photons/plot_isr_photons.sh
```

The existing recipe reads only the five Florian `photosFSR` samples (`Zee`,
`Zmumu`, `Ztautau`, `ZKK`, `Zpipi`), totaling 449,999 events. `python/delphi_analysis/photons/data.py`
resolves truth, and `plot_isr_photons.py` produces the same plot set for all groups:

```text
plots/20260828_florian/photon_study/
  01_gamma/<sample>/       all stable Gen photons (PDG 22, status 1)
  02_ISR/<sample>/         stable ISR, including Beam ISR
  03_noBeamISR/<sample>/   stable ISR excluding Beam ISR
  study_summary.json      event, topology, association, and efficiency counts
```

The existing ISR ancestry definition is unchanged: the first non-photon ancestor
is an electron/positron, with no hard parent on the path. Beam ISR additionally
has exactly zero px and py and a direct, parentless incoming electron/positron
parent. No energy, fiducial-angle, or angular-matching cut is applied.
Previous results remain under `plots/20260828_florian/isr_photons/`.

**One Gen photon is the analysis unit.** For every saved SimPart, follow
`SimPart_genIdx`, or, when negative, trace parents through
`SimPart_originVtxIdx → SimVtx_incomingIdx` to the first valid Gen anchor.
For reco candidates, enter this graph via `Photon_partIdx → Part_simIdx`.
Any Sim species can carry the reco association; several reco fragments count
as just one reconstructed Gen photon. Angular matching is neither run nor used.

Topology is operational and mutually exclusive:

- **no_sim:** no identifiable saved Sim lineage.
- **direct:** saved lineage, but no secondary Sim descendants.
- **shower:** at least one saved secondary Sim descendant, of any species.

Here `shower` means stored secondary activity; it does not identify a particular
conversion/shower process or prove a detector energy deposit. Sim γ multiplicity
counts **secondary photons only**, excluding directly Gen-linked anchors.
`SimPart_pdgId` contains DELPHI mass codes: photons are code **21**, not PDG 22
([SKELANA manual, printed page 42](https://opendata.cern/record/80502/files/skelana.pdf#page=44)).
A conversion with only electron descendants is still shower topology even when
its secondary Sim γ count is zero.

Let G be selected Gen photons, S those with saved Sim lineage, and R those with
at least one associated reco Photon:

```text
Gen → Sim = N(S)/N(G)
Sim → Reco = N(R)/N(S)   conditional, still counted per Gen photon
Gen → Reco = N(R)/N(G)
```

All three use the same **Gen E and Gen cos(theta)** bins, with 68.27%
Clopper–Pearson intervals. The code checks topology coverage, `R ⊆ S ⊆ G`,
and binwise central-value closure `εGR = εGS × εSR`. Empty denominators are
undefined, not zero. The three ratios are correlated.

Every group has Gen energy/cos(theta)/2D distributions, topology fractions,
separate GS/SR/GR efficiencies, Sim-γ/reco multiplicities and their 2D correlation,
and leading/summed energy ratios and leading-reco opening angle. Response
uses only matched Gen photons; leading means highest reco energy. Multiplicity
is normalized to all Gen photons, response to matched Gen photons, and Gen
spectra to events. Count plots have sqrt(N) errors. All ranges retain all entries.
Logarithmic figures have fully linear `_linear.png` siblings; cos(theta)
fraction/efficiency figures are already linear. The angular-response symlog
transition at 0.001 rad is only a display scale. CMS style, default font sizes,
11×9 figures, DELPHI Simulation, and bold sample legends are retained.

In these samples every stable Gen photon has a saved Sim anchor, so GS=100%
and SR=GR. This is **saved-lineage coverage, not detector transport efficiency**.
Beam ISR has no stored secondary descendants or associated reco candidates.
Candidates without stored truth remain unassociated; the summary reports them.
The saved schema cannot separate transport, acceptance, or missing-association
causes, and no detector-interaction position is inferred from `SimVtx_position`.

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
