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
isr-photons.py                  Gen-photon reconstruction and associated energy response
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
  study_summary.json      efficiency, energy response, and truth-link diagnostics
  manual_audit/           preserved event trees and manual review
```

The existing ISR ancestry definition is unchanged: the first non-photon ancestor
is an electron/positron, with no hard parent on the path. Beam ISR additionally
has exactly zero px and py and a direct, parentless incoming electron/positron
parent. No energy, fiducial-angle, or angular-matching cut is applied.
Previous results remain under `plots/20260828_florian/isr_photons/`.

**One stable Gen photon is the analysis unit.** Associate each Reco Photon via
`Photon_partIdx → Part_simIdx → SimPart_genIdx`; when the latter is negative,
trace Sim parents through `SimPart_originVtxIdx → SimVtx_incomingIdx`.
Stop at the **first Gen anchor**, even when it has an earlier Gen parent.
Associate the candidate only if that origin is a selected stable Gen photon.
Direct and descendant Sim links are equally valid; the linked Sim need not be
a photon or a terminal node. There is no angular fallback. Candidates without
usable truth remain unresolved; valid origins outside stable photons are counted
separately in the summary.

```text
Gen → Reco = N(Gen photons with >=1 associated Reco Photon) / N(selected Gen photons)
Energy recovery per Gen photon = sum(associated Reco Photon energies) / E_gen
```

Several reco candidates count as **one efficiency success**, while all their
energies enter the sum. Unmatched photons have recovery **zero**, not an undefined
value. Leading response uses the highest-energy candidate and is conditional on
successful reconstruction. Recovery ratios are not clipped at one.

Each sample/population has the same nine observables:

| Plots | Meaning |
| --- | --- |
| `gen_to_reco_efficiency_vs_{energy,cos_theta}` | Gen-unit success fraction, 68.27% Clopper–Pearson intervals |
| `summed_energy_response` | Recovery distribution over all Gen photons, including unmatched zeros |
| `energy_recovery_vs_{energy,cos_theta}` | Unweighted mean recovery per Gen photon, including zeros; standard error of the mean |
| `leading_energy_response` | Highest-energy reco response among matched Gen photons |
| `reco_gamma_multiplicity` | Reco candidates per Gen photon, including zero |
| `linked_sim_depth` | Parent steps from the reco-linked Sim to its first Gen anchor; zero is direct |
| `linked_sim_species` | Linked Sim species, stacked as direct/descendant |

The last two are normalized per associated **Reco Photon**; the other distributions
are normalized per Gen photon (leading response per matched Gen photon). Count
errors are sqrt(N). Empty efficiency/profile bins are undefined. A singleton
profile bin has no estimable standard error. Energy plots share binning across
samples within each population: about three bins/decade, or 1 GeV for linear
figures; cos(theta) uses 20 bins. Histograms retain all entries and tails.
Logarithmic figures retain fully linear `_linear.png` siblings (16 files per
sample/population). Summed response uses a linear interval at zero followed by
a logarithmic x-axis; its first bin contains unmatched zeros only. Linear response
histograms span the full tail and may merge the central response into a broad bin.
CMS style, default font sizes, 11×9 figures, DELPHI Simulation, and bold sample
legends are retained.

Sim coverage and saved-secondary counts remain summary diagnostics. Separate
GS/SR efficiencies, topology fractions, Sim-gamma multiplicity/correlation,
standalone Gen spectra/maps, and angular response plots are removed. Every stable
photon in these samples has a saved Sim anchor; this is **saved-lineage coverage,
not detector transport efficiency**. `SimPart_pdgId` stores DELPHI mass codes
(gamma=21, electron/positron=±2), not PDG IDs. Beam ISR has no associated reco.

**Energy attribution is a limitation of this truth record.** A single `Part_simIdx`
does not provide the energy contribution of every Sim particle to a reco candidate.
Very soft Gen photons can be linked to much more energetic candidates (audited ISR:
0.000891 GeV Gen → 1.100 GeV reco). The observable is therefore **associated reco
energy response**, not proof that this energy came exclusively from that photon.
Large tails can dominate arithmetic means; matched 16/50/84% quantiles and the
ratio of total associated reco to total Gen energy are also saved in the summary.
No energy-consistency cut is silently applied. Missing truth links and incomplete
saved trees prevent separating transport, acceptance, and association failures:
362,929 of 573,796 reco candidates have no usable Part→Sim link. Thus the efficiency
measures saved truth-associated reconstruction, rather than proving all physically
reconstructed photons have been identified. `study_notes.txt` records validation
results and the two extreme-response event examples.

For an event-level truth audit, run:

```bash
bash runs/photons/dump_truth_trees.sh
```

This saves ten event dumps (two per sample) and a selection index in
`plots/20260828_florian/photon_study/manual_audit/`. Each dump includes the full
Gen record, the full stored Sim tree, resolved Gen origins, and every reco Photon
association. A baseline and a targeted topology case are selected per sample;
there are no energy or angular cuts. The script checks parent vertices, tree
coverage, Gen anchors, and agreement between forward and backward lineage tracing.
`REVIEW.txt` in that output folder records the manual review of the selected events.

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
