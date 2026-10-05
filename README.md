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
`Zpipi`) are read, totaling 449,999 events. Code and recipes live under
`python/delphi_analysis/isr_study/`, `scripts/isr_study/`, and `runs/isr_study/`.
The plotting CLI has only input/output directory arguments.

```text
plots/20260828_florian/isr_study/
  00_link_validation/
    <sample>/                        link agreement and conflicts
    manual_audit/                    representative full event trees
  01_photon_conversion/
    gen/<sample>/                    stacked Beam/non-Beam Gen ISR spectra
    {all_ISR,noBeamISR}/<sample>/     Photon or reconstructed conversion
  02_all_lineage/
    {all_ISR,noBeamISR}/<sample>/     all associated reconstructed Part objects
  03_angle_matching/
    01_photon_conversion/             angular cross-check for photon-like objects
    02_all_lineage/                    angular cross-check for all Part activity
  study_summary.json                 counts, response summaries, and diagnostics
  study_notes.txt                    assumptions and unresolved limitations
```

**One stable Gen ISR photon is the analysis unit:** `GenPart_pdgId == 22`,
`GenPart_status == 1`, first non-photon Gen ancestor an electron/positron, and no
hard-process parent along that path. Beam ISR additionally has exactly zero px
and py and a direct, parentless incoming electron/positron parent. Results are
separated into `all_ISR` and `noBeamISR`; Gen spectra show Beam and non-Beam ISR
as stacked components. No energy, fiducial-angle, PID, lock, or angular cut is
applied to nominal truth reconstruction.

### Stored truth association

Compare both saved directions independently:

- `Part_simIdx` -> Sim ancestry -> first existing `SimPart_genIdx` anchor;
- directly Gen-linked Sim anchors -> descendants -> `SimPart_partIdx`.

A nested Gen anchor starts a different lineage. Direct and descendant Sim objects
of any species are valid; the linked object need not be a gamma or a terminal
node. All known Gen origins, including non-ISR origins, are retained when checking
agreement. A Part is accepted if the **union contains exactly one Gen origin**.
Missing evidence in one direction is allowed; contradictory or absent origins
remain unresolved. Angular matching does not repair these nominal links.

### Photon/conversion and all-lineage channels

`01_photon_conversion` measures the fraction of Gen ISR photons with at least one
associated Photon candidate or reconstructed conversion. `02_all_lineage`
measures the fraction with at least one associated Reco Part. Several candidates
still count as one Gen success. Efficiencies and mean energy responses are shown
versus Gen energy, pT, and cos(theta), for both populations.

The `Photon` branch is a view of neutral Parts with EM calorimeter energy, not a
pure photon PID selection. `Part_pdgId` and `SimPart_pdgId` contain DELPHI mass
codes, not PDG IDs (gamma=21, electron/positron=±2). Full-lineage composition uses
coarse reconstructed charge/mass-code categories; these are not perfect truth PID.
Composition and multiplicity count raw associated Parts, including retained parent
and daughter representations.

`PhotonConv_simPhotonIdx` is unusable in these samples (all values are -1).
Conversions instead use their saved daughter origin vertex and its incoming Part.
A neutral gamma-code parent with a reciprocal decay vertex validates the topology.
The conversion is accepted only when all known parent/daughter Gen anchors agree.
One available anchor can suffice; a known conflicting anchor vetoes association.

For energy response, Photon rows and their source Parts are the same representation.
A validated conversion represents its parent and daughters once; it is not summed
again with a Photon view of that group. Full-lineage response keeps an associated
composite parent instead of its descendants when their known origins are consistent.
Duplicate conversion rows sharing a saved parent count once; an accepted ancestor
conversion also replaces a nested conversion of the same consistent lineage.
Parts that merely share a Sim link are not deduplicated without saved
parent/daughter structure.
Accepted objects in groups with contradictory anchors have undefined energy
response; raw Part activity and multiplicity remain separate diagnostics. Undefined
response is not replaced by zero. A Gen photon with no accepted channel candidate
has zero associated response, including rejected conversion-only evidence.
Ratios are not clipped at one.

Plots retain efficiency, summed/leading response, response profiles, and multiplicity.
Full-lineage plots additionally show raw activity response and stacked composition
versus pT and cos(theta). Response profiles use defined responses only; numerical
summaries record undefined-response coverage.

### Angular cross-check and limitations

`03_angle_matching` scans 3D opening-angle cuts separately for the two channels.
Truth lineage supplies the reference: own-ISR recovery, known unrelated origins,
and candidates with unresolved truth are distinguished. Unknown truth is not
labelled a false positive. Energy recovery and candidate multiplicity are retained
in the scan; no preferred cut is silently promoted to nominal matching. The cone
check allows several candidates per Gen photon; it is not an exclusive assignment.
The recovery curve is conditional on nominal truth success; the any-candidate
curve uses all selected Gen ISR photons. Pair fractions count Gen–Reco pairs,
including repeated candidates in overlapping cones.
The cone scan uses the deduplicated candidate representations. Full-lineage
nominal success and composition still count the raw Part footprint.
Energy profiles omit cones containing ambiguous representations, with coverage
recorded in the summary. An unresolved saved parent retained alongside its daughters
is energy-invalid in this geometric check; nominal daughter-associated energy and
truth success are unchanged. Zero-momentum candidates are excluded only from
geometric cones because their direction is undefined, and are counted in diagnostics.

The efficiencies measure **saved truth-associated reconstruction**, not detector
transport efficiency. Saved Gen-to-Sim coverage is a link diagnostic. Missing
links, retained rejected Parts, and mixed contributions to one reconstructed object
remain limitations. A single origin label does not measure its fractional energy
contribution: very soft ISR can be associated with much larger reco energy. Even
a deduplicated response is therefore an **associated energy response**, not proof
that all that energy originated from the selected photon.

Count errors use sqrt(N); efficiencies use 68.27% Clopper–Pearson intervals.
Empty efficiency/profile bins are undefined. Profiles use the standard error of
the mean where estimable. CMS style, default font sizes, 11×9 figures, DELPHI
Simulation, and bold sample legends are retained. Logarithmic figures have fully
linear `_linear.png` siblings where useful.

The truth-tree script preserves the ten original representative events and adds
fixed conversion, reverse-only, contradictory-anchor, and suspicious Beam ISR
cases. Dumps include the full Gen/Sim records, Photon comparison, all Part link
origins, and conversion associations. Independent checks cover Sim-tree structure,
nested-anchor boundaries, and the conflict-free Part union. Output is under
`00_link_validation/manual_audit/`; previous study outputs are preserved under
`plots/20260828_florian/archive/`.

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
