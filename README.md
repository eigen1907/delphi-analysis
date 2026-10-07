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
  isr_study/          ISR links, conversions, full reco lineage, and angle checks
scripts/{checks,tracking,rich,pid,isr_study,data}/
runs/{checks,tracking,rich,pid,isr_study,data}/
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
bash runs/isr_study/plot_isr_study.sh
bash runs/pid/train_bdt_pid_standard.sh
bash runs/pid/apply_bdt_pid_standard.sh
```

The data preparation, checks, and BDT recipes target `20260606_100kTest`
under `data/202606xx_jongwon/`. The ISR recipe targets the Florian samples.
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
isr_study/plot.py               ISR reconstruction, conversions, and angular cross-checks
```

Plotting diagnostics and CSV summaries are written below
`data/check/<sample-set>/`.

## ISR reconstruction study

```bash
bash runs/isr_study/plot_isr_study.sh
```

The five Florian `photosFSR` samples (`Zee`, `Zmumu`, `Ztautau`, `ZKK`, `Zpipi`)
contain 449,999 events. The plotting script accepts only input/output directories.
Code and recipes are grouped under `isr_study/`.

```text
plots/20260828_florian/isr_study/
  {gamma,gamma_plus_conversion,gamma_pid}/
    01_sample_distribution/{gen_gamma,gen_isr,gen_no_isr,gen_isr_non_beam}/
    02_truth_link_matching_validation/
    03_angular_matching_validation/
    04_matching_result/{truth_matching,angular_matching}/{no_cut,0p1_cut}/
    05_detector_efficiency/{truth_matching,angular_matching}/{no_cut,0p1_cut}/
```

Stage 01 keeps four stable Gen selections, including Beam ISR for inspection:

| Selection | Definition |
| --- | --- |
| `gen_gamma` | All stable gamma |
| `gen_isr` | All ISR |
| `gen_no_isr` | Gamma excluding ISR: FSR + Decayed |
| `gen_isr_non_beam` | ISR excluding Beam ISR |

Stable photons have `GenPart_pdgId == 22` and `GenPart_status == 1`.
ISR has an electron/positron as its first non-photon ancestor, with no hard parent
on that path. Beam ISR additionally has `px == py == 0` and a direct, parentless
incoming electron/positron parent. FSR starts at the hard parent or its direct
sample-final-state daughter. Decayed is the previous Others category; here it
contains π⁰, η, and ω decay photons in Ztautau. Stages 04 and 05 compare the full
Gen energy range (`no_cut`) with Gen E ≥ 0.1 GeV (`0p1_cut`) for every plot.
No additional angular, fiducial, or lock cuts are applied.

Stages 02–05 exclude Beam ISR from all Gen selections, numerators, and
denominators, and compare three populations in each plot:

| Legend | Definition |
| --- | --- |
| $\gamma_{\mathrm{all}}$ | Stable gamma excluding Beam ISR |
| $\gamma_{\mathrm{ISR}}$ | Non-beam ISR |
| $\gamma_{\mathrm{no\ ISR}}$ | FSR + Decayed |

This exclusion uses the Beam ISR ancestry and exact `px == py == 0` definition;
it is not an energy threshold or detector acceptance cut. It changes only the
study populations: the original Gen particles and truth ancestry remain intact.

The same study runs for three Reco collections: **Photon only** (`gamma`),
**Photon + conversion** (`gamma_plus_conversion`), and **Part gamma mass code**
(`gamma_pid`). Photon only retains every original Photon row. The combined
collection replaces overlapping Photon representations once by the conversion
candidate. `gamma_pid` selects **all Parts with `Part_pdgId == 21`**, uses
`Part_fourMomentum`, and appears as gamma_Part in Reco legends. It does not
require membership in Photon or add PhotonConv candidates. All three use the
same Gen populations, bins, and truth definitions, reading each event once. The Photon
branch is a view of neutral Parts with EM calorimeter energy, not pure photon PID.
`Part_pdgId` and `SimPart_pdgId` use DELPHI mass codes: gamma=21, e±=±2.
The `gamma_pid` name describes the stored Part classification, not a strict
photon-ID working point: the DST writer also assigns code 21 to unidentified
neutral Parts. The original identification before that fallback is not saved
separately. `Part_massId` is a packed identification/topology word, so requiring
`massId == 21` or `massId != 0` would not recover an explicit gamma PID. No such
condition, charge cut, or lock cut is added. Code-21 Parts can include conversion
parents and neutral objects outside the Photon view. Overlapping saved
parent/descendant representations retain their counts but have undefined energy
sums, as in the Photon-only study.

- **01:** Stacked Gen colors and one inclusive Reco outline with
  hatching on one axis. Reco is identical across the four Gen selections. Energy and
  cos(theta) have `_log` versions. Linear energy uses a labelled E ≥ 50 GeV
  overflow bin; logarithmic energy retains the full tail. Each Gen directory has
  event-multiplicity plots with its own Gen components: all four in `gen_gamma`,
  BeamISR/NonBeamISR in `gen_isr`, FSR/Decayed in `gen_no_isr`, and NonBeamISR in
  `gen_isr_non_beam`. One inclusive Reco curve appears in every plot.
  Multiplicity axes are N_gamma and raw Events;
  `_log` changes only the multiplicity y-axis.
- **02:** One plot per sample compares the three Gen selections, with Agree,
  Forward only, Reverse only, Conflict, and No origin on the x-axis. Only candidates
  with known evidence involving that selection enter the counts; No origin cannot
  be assigned to a Gen selection and is therefore zero. Status is evaluated using
  all original anchors before applying the Gen selection, so a conflict between
  Beam ISR and another origin remains a Conflict. Conflict is a link diagnostic,
  not an accepted match.
- **03:** The three Gen selections scan the 3D opening angle, displayed to 10°.
  This truth-retention scan uses the full Gen energy range. Curves retain each
  Gen photon's own truth-associated candidates. For each Reco
  collection, the reference cut retains ≥99% of pooled non-beam ISR truth successes;
  it is a truth-retention check, not a purity optimum or a nominal truth cut.
- **04:** Both `truth_matching` and `angular_matching` contain `no_cut` and
  `0p1_cut` subdirectories, each saving the same six plots per sample. `no_cut`
  has no Gen energy threshold; `0p1_cut` requires **Gen E ≥ 0.1 GeV** in every
  plot, including both the efficiency numerator and denominator. Beam ISR remains
  excluded in both versions. Each plot compares All, ISR, and no ISR. The plot set
  is efficiency versus Gen E and cos(theta), plus two energy-ratio histograms:

  - **Per Gen photon:** `sum(E_reco) / E_gen`, requiring at least one associated
    Reco candidate and a finite, unambiguous energy sum. One matched Gen photon
    contributes one entry; unmatched photons remain in the efficiency denominator.
  - **Per event:** `sum(E_reco) / sum(E_gen)`. The denominator includes every
    selected Gen photon, including unmatched ones. The numerator sums distinct
    associated Reco candidates once; angular matching uses the union of the
    selected Gen cones, without a truth gate. Events with selected Gen photons
    but no associated Reco enter at zero. Events with no selected Gen photons or
    an ambiguous energy sum are excluded.

  Efficiency counts ≥1 candidate once per Gen photon, even when its energy sum
  is ambiguous. Both ratio histograms span 0–3 in 0.1 bins, without percent
  scaling, with count errors and linear-y/`_log` versions. Ratio one means equal
  energy. Their filenames are `matched_energy_ratio_<sample>[_log].png` and
  `event_energy_ratio_<sample>[_log].png`.
  Separate one-column legends align Gen selection and statistics, using empty
  handles. The matched-photon mean uses unbinned ratios **0 ≤ R ≤ 3**; the event
  mean uses **0 < R < 3**, excluding both endpoints. Histograms still include
  zero and three. The event legend also reports `Frac. (=0)`. Zero and >3 tail
  fractions use the same denominator: all eligible finite event ratios, including
  zero and overflow. Photon tail fractions use all finite matched-photon ratios.
  Tails are neither folded into bins nor clipped; excluded and tail counts are
  printed to the terminal. No Reco energy cut is added.
- **05:** Five-sample pooled regional efficiency uses the same two matching
  methods, two Gen energy versions, and three Gen populations for each Reco mode.
  Each case has a separate plot with counts and Clopper–Pearson intervals. The
  Gen threshold applies to both numerator and denominator; Beam ISR is excluded.
  All angles is
  compared with HPC (40–140°), FEMC (10–37° and 143–170°), and STIC (2–10° and
  170–178°), using open angular intervals. No Reco energy cut is applied.
  A Gen photon succeeds once when it has at least one matched Reco candidate,
  even if its associated energy sum is ambiguous. Truth results use the saved
  association; angular results use the stage 03 cone cut with no truth gate.
  Missing forward truth links and accidental angular associations limit a direct
  interpretation as detector performance. Filenames are
  `efficiency_by_detector_<population>_combined.png`.

Nominal truth compares `Part_simIdx → Sim ancestry → first Gen anchor` with
`Gen-linked Sim → descendants → SimPart_partIdx`. A unique union origin is
accepted; absent or contradictory origins remain unresolved. A nested Gen anchor
starts a new lineage. Direct and descendant Sim objects of any species are valid;
angular matching never repairs truth. `PhotonConv_simPhotonIdx` is always -1 here,
so conversions use reciprocal saved parent/daughter vertices and agreeing anchors.

Angular cones can share candidates and are not exclusive assignments. Zero-momentum
Reco objects have no angular direction. Overlapping saved representations invalidate
an energy sum, including parent/descendant Photon rows in the Photon-only study;
candidate counts and efficiency still retain those rows. Responses can exceed one
because candidates may contain mixed energy, especially for soft Gen photons.
The per-photon ratio histogram is conditional on a match; the event ratio also
includes energy from unmatched Gen photons in its denominator. Angular event sums
count each Reco candidate once, even when it enters several Gen cones. A truth
link identifies ancestry, not exclusive ownership of the Reco energy.
Saved Gen-to-Sim coverage is a link diagnostic,
not detector transport efficiency; missing links and reverse-only Beam associations
remain limitations. No summary JSON, note, or validation text files are produced.

CMS/mplhep style and DELPHI Simulation are retained, with slightly thicker curves.
Legends use 22 pt text in one column at the upper right without a frame; other
text keeps the CMS defaults. A bold symbolic sample label appears separately
at the upper left, e.g. Z → μ⁺μ⁻. Only `0p1_cut` plots display the Gen energy cut
caption. No detector-region guides are drawn.
Counts use sqrt(N) errors; efficiencies use 68.27% Clopper–Pearson intervals.
Linear energy bins start at 1 GeV width and widen at higher energies; this adds
no energy selection. Log-energy plots retain equally spaced logarithmic bins.
Filenames place the sample before the scale suffix: `energy_Zee_log.png`.

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
