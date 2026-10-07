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
    04_matching_result/{truth_matching,angular_matching}/
    05_detector_efficiency/
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
contains π⁰, η, and ω decay photons in Ztautau. Efficiency curves retain the full
Gen energy range. Stage 04 energy-response diagnostics and the stage 05 regional
summary use Gen E ≥ 0.1 GeV. No additional angular, fiducial, or lock cuts are applied.

Stages 02–04 exclude Beam ISR from all Gen selections, numerators, and
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
- **04:** Separate `truth_matching` and `angular_matching` directories contain
  efficiency and energy response versus Gen E and cos(theta). Each plot compares
  the three Gen selections. Angular matching retains the inclusive Reco pool,
  without a Reco truth-origin veto, and includes every candidate inside the
  chosen cone. Efficiency counts ≥1 candidate once per Gen photon. Response is
  sum(E_reco)/E_gen, including unmatched zeros; ambiguous sums are undefined.
  Response profiles require **Gen E ≥ 0.1 GeV**; the efficiency curves have no
  energy threshold. Response plots also have `_log` versions.
  No matched-multiplicity plots remain.
  Matched energy distributions use one Gen photon as the unit:
  `(E_gen - sum(E_reco)) / E_gen` in **−2 to +2**, and
  `sum(E_reco) / E_gen` in **0 to 3**, without percent scaling. Residual zero
  and ratio one both mean equal energy; positive residual means energy loss,
  and negative residual means excess Reco energy. Both distributions require
  **Gen E ≥ 0.1 GeV**, at least one associated Reco candidate,
  and a finite, unambiguous energy sum. Unmatched photons remain in the efficiency
  denominator but do not enter this matched-response diagnostic. Bins have width
  0.1 and count errors. Reco energies are nonnegative, so residuals cannot exceed
  one: the +1 to +2 interval is intentionally empty. The two plots show the same
  information, related by residual = 1 − ratio.
  Separate one-column legends align Gen selection, `Mean (shown)`, and the tail
  fraction, using empty handles for statistics. The mean uses unbinned values
  inside the displayed range, including its endpoints. Residual < −2 and ratio > 3
  are reported as decimal fractions of all finite matched photons passing the
  energy cut. Out-of-range values are neither folded into bins nor clipped.
  Counts of ambiguous energy sums and tails are printed to the terminal.
  Linear-y and `_log` versions are saved as
  `matched_energy_{residual,ratio}_<sample>[_log].png`.
  Gen-versus-Reco energy histograms use the same finite matched selection:
  x = E_gen and y = sum(E_reco), one entry per Gen photon. All, ISR, and no ISR
  have separate images within each matching directory, named
  `gen_reco_energy_2d_<population>_<sample>[_log].png`. Both energy axes are linear;
  `_log` uses a logarithmic count color scale. Empty cells are blank and the
  dashed diagonal denotes equal energies. All Reco sums ≥ 50 GeV enter the
  labelled overflow bin, so large responses are retained rather than dropped.
  The first Gen energy bin starts at 0.1 GeV and ends at 1 GeV; the threshold
  selects Gen photons, without adding a Reco energy requirement.
- **05:** Five-sample pooled **truth-matched** non-beam ISR efficiency for
  **Gen E ≥ 0.1 GeV**, with counts and Clopper–Pearson intervals. All angles is
  compared with HPC (40–140°), FEMC (10–37° and 143–170°), and STIC (2–10° and
  170–178°), using open angular intervals. No Reco energy cut is applied. This measures
  reconstruction with usable saved truth association. Missing truth links among
  forward Reco candidates limit the interpretation of the very low STIC fraction
  as a measure of physical detector performance.

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
Response profiles average per-Gen ratios, not total recovered energy divided by
total Gen energy. The stage 04 Gen energy cut limits sensitivity to very soft photons with large ratios;
a truth link identifies ancestry, not exclusive ownership of the Reco energy.
Saved Gen-to-Sim coverage is a link diagnostic,
not detector transport efficiency; missing links and reverse-only Beam associations
remain limitations. No summary JSON, note, or validation text files are produced.

CMS/mplhep defaults and DELPHI Simulation are retained. Legends use one column
at the upper right without a frame; a bold symbolic sample label appears separately
at the upper left, e.g. Z → μ⁺μ⁻. No detector-region guides are drawn.
Counts use sqrt(N) errors; efficiencies use 68.27% Clopper–Pearson intervals.
Linear energy bins start at 1 GeV width and widen at higher energies; this adds
no energy selection. Log-energy plots retain equally spaced logarithmic bins.
Response profiles show the standard error of the mean where estimable. Filenames
place the sample before the scale suffix: `energy_Zee_log.png`.

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
