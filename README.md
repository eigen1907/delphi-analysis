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
bash runs/isr_study/dump_truth_trees.sh
```

Only the five Florian `photosFSR` samples (`Zee`, `Zmumu`, `Ztautau`, `ZKK`,
`Zpipi`) are read: 449,999 events. Code and recipes remain under
`python/delphi_analysis/isr_study/`, `scripts/isr_study/`, and `runs/isr_study/`.
The plotting CLI has only input/output directory arguments.

```text
plots/20260828_florian/isr_study/
  01_sample_distribution/
    {gen_gamma_all,gen_isr_all,gen_isr_non_beam}/
    multiplicity/
  02_truth_link_matching_validation/
    <gen_selection>/<reco_definition>/
    manual_audit/
  03_angular_matching_validation/
  04_matching_result/
    {gen_all_gamma,gen_all_isr,gen_non_beam_isr}/
  study_summary.json
  study_notes.txt
```

All figure filenames end in the sample name. CMS style, default font sizes,
DELPHI Simulation, and bold sample legends are retained. Count errors use
sqrt(N); efficiencies use 68.27% Clopper–Pearson intervals. Empty bins are
undefined; response profiles show the standard error of the mean where estimable.
Energy plots use wider bins at high energy, with a separate logarithmic-energy
version. Linear sample spectra use a labelled overflow bin for E >= 50 GeV;
the log-energy plots and matching use the original energies and full tails.
No pT plots are produced. Response and matched-multiplicity plots have
additional `_logy` versions.

### Gen populations and Reco definitions

**One stable Gen photon is the analysis unit:** `GenPart_pdgId == 22` and
`GenPart_status == 1`. The three selections are all stable photons, all ISR,
and non-beam ISR. ISR has an electron/positron as its first non-photon Gen
ancestor, with no hard-process parent on that path. Beam ISR additionally has
exactly zero px and py and a direct, parentless incoming electron/positron parent.
FSR has the hard-process parent or its direct sample-final-state daughter as the
first non-photon ancestor. `Decayed` is the previous `Others` category; in these
samples its 197,201 photons come exclusively from π⁰, η, and ω decays in Ztautau.
No energy, fiducial-angle, PID, lock, or angular cut is applied to nominal truth
reconstruction.

Reco definitions are Photon candidates (`reco_gamma`), Photon candidates plus
explicit reconstructed conversions (`reco_gamma_plus_conv`), and all associated
Parts (`reco_all`). The `Photon` branch is a view of neutral Parts with EM
calorimeter energy, not a pure photon PID selection. `Part_pdgId` and
`SimPart_pdgId` use DELPHI mass codes: gamma=21, electron/positron=±2.

`01_sample_distribution` shows Gen stacks above separate Reco gamma/conversion
stacks, versus energy and cos(theta). Gen stacks retain all selected photons;
Reco ISR subsets require a unique selected truth origin. Reco candidates with
unresolved origins cannot be assigned to ISR subsets. The three Gen multiplicity
populations count matched raw Parts **per Gen
photon**, including zero matches. The two Reco multiplicity populations count
gamma or gamma-plus-conversion candidates **per event**. The multiplicity in
`04_matching_result` compares associated Reco objects **per Gen photon** for
all three Reco definitions.

### Stored truth association

Compare both saved directions independently:

- `Part_simIdx` → Sim ancestry → first existing `SimPart_genIdx` anchor;
- directly Gen-linked Sim anchors → descendants → `SimPart_partIdx`.

A nested Gen anchor starts a different lineage. Direct and descendant Sim objects
of any species are valid; the linked object need not be a gamma or a terminal
node. All known Gen origins, including non-photon origins, are retained when
checking agreement. A Part is accepted if the **union contains exactly one Gen
origin**. Missing evidence in one direction is allowed; contradictory or absent
origins remain unresolved. Angular matching never repairs nominal truth links.

`02_truth_link_matching_validation` contains the nine Gen-selection × Reco-definition
combinations, with Agree, Forward only, Reverse only, Conflict, and No origin.
Each plot counts Reco objects with known evidence involving the selected Gen
population, including conflicts. Its No origin bin is therefore zero by
construction. The summary separately retains all-Reco unresolved counts and
known origins outside each Gen selection; unknown origins are not invented to
fill a selected population.
These categories compare saved origin evidence; conversion acceptance additionally
requires the saved parent/daughter topology described below.

`PhotonConv_simPhotonIdx` is unusable here (all values are -1). Conversions use
their daughter origin vertex and its incoming Part. A neutral gamma-code parent
with a reciprocal decay vertex validates the topology. Association requires all
known parent/daughter Gen anchors to agree; one available anchor can suffice.
Duplicate conversion rows sharing a parent count once, as do consistent nested
conversion representations.

### Efficiency, energy response, and angular validation

`04_matching_result` compares the three Reco definitions in each Gen population.
Efficiency is the fraction of Gen photons with at least one associated Reco
object; several objects still count as one success. Energy response is
`sum(E_reco) / E_gen`, shown versus Gen energy and cos(theta).

Photon rows and their source Parts are the same representation. A validated
conversion represents its parent and daughters once, without adding a Photon
view of the same group. For all-Part energy, an associated composite parent
replaces consistent saved descendants. Merely sharing a Sim link is insufficient
for deduplication. All-Part multiplicity still counts raw associated Parts.
Ambiguous associated energy is NaN, never zero; no accepted candidate gives zero
response. Profiles use defined responses, including unmatched zeros, and the
summary records coverage. Ratios are not clipped at one. These are **associated
energy responses**, not measurements of the selected photon's pure energy
contribution: reconstructed objects can contain mixed contributions, especially
for very soft Gen photons.

`03_angular_matching_validation` compares six Gen/Reco combinations per sample:
the three Gen populations with gamma or gamma plus conversion. A 3D opening-angle
cone can contain several candidates and overlapping cones can reuse candidates.
The main curves count only each Gen photon's own truth-associated candidates
inside the cone: efficiency uses all selected Gen photons as denominator, and
energy recovery sums their cone energy divided by Gen energy, including
unmatched zeros. All-candidate cone results, known unrelated origins, and
unresolved candidates remain diagnostics in the JSON summary. Unknown truth is
not automatically labelled a false positive.
The dashed cut is the smallest scanned angle retaining at least 99% of nominal
truth successes for pooled non-beam ISR in **both** photon channels. It is a
truth-retention cross-check, not a purity optimum or a cut on nominal results.
Zero-momentum Reco candidates are excluded only from angular cones because their
direction is undefined. Cones with structurally ambiguous energy are omitted
from energy profiles, with coverage recorded in the summary.

The efficiencies measure saved truth-associated reconstruction, not detector
transport efficiency. Missing links, contradictory anchors, and suspicious
reverse-only or soft-Beam associations remain limitations. Saved Gen-to-Sim
coverage is a link diagnostic rather than a main physics efficiency.

The manual audit preserves ten original representative trees plus eight fixed
conversion, reverse-only, contradictory-anchor, and suspicious Beam ISR cases.
Dumps include full Gen/Sim records, Photon comparisons, all Part origins, and
conversion associations. Independent checks cover Sim-tree structure,
nested-anchor boundaries, and the conflict-free Part union. Output is under
`02_truth_link_matching_validation/manual_audit/`; previous outputs are preserved
under `plots/20260828_florian/archive/`.

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
