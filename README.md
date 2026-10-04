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
`Ztautau`, `ZKK`, and `Zpipi` (449,999 events). Implementation is in
`python/delphi_analysis/photons/`: `data.py` resolves origins,
`plot_stages.py` draws the Gen–Sim–Reco decomposition,
`plot_isr_photons.py` runs the study and geometric diagnostics, and
`plot_response.py`, `plot_2d.py`, and `plot_generator.py` contain ISR checks.

### Three photon groups

Plots are saved below `plots/20260828_florian/isr_photons/`:

```text
01_gen_sim_reco/<sample>/              all stable gen photons
02_gen_sim_reco_ISR/<sample>/          all stable gen ISR, including Beam ISR
03_gen_sim_reco_nobeamISR/<sample>/    stable gen ISR excluding Beam ISR
04_geometric_matching/<sample>/      auxiliary angular-matching comparisons
04_geometric_matching/                pooled angle-cut validation
study_summary.json                    counts, efficiencies, and bookkeeping checks
```

The previous plots are preserved in the sibling directory
`isr_photons_previous_20261004/`; the recipe produces the new layout.

Stable gen photons have `GenPart_pdgId == 22` and `GenPart_status == 1`.
The existing gen ancestry rules are unchanged: the first non-photon ancestor
is an electron/positron and the path contains no hard parent for ISR; the hard
parent or its sample-specific final-state daughter identifies FSR. Other
origins are Others. Beam ISR is ISR with exactly zero px and py and a direct,
parentless incoming electron/positron parent. Non-beam ISR is the remainder.
There is **no energy threshold, fiducial-angle cut, or angular/energy matching
cut** in the nominal study. Invalid photon four-vectors stop the run rather
than being removed. Undefined phi is displayed at zero, including Beam ISR.

### Gen–Sim–Reco efficiencies

Each group uses the same unit: one unique selected **gen photon**.

- **G:** all selected gen photons.
- **S:** gen photons with at least one saved Sim object resolving to that Gen origin.
- **R:** gen photons with at least one reco Photon associated with that Sim lineage.

For every Sim object, use `SimPart_genIdx` when nonnegative; otherwise follow
`SimPart_originVtxIdx → SimVtx_incomingIdx` until the first valid Gen index.
For reco, use `Photon_partIdx → Part_simIdx` to enter the same Sim graph.
The direct Sim root and descendants of any particle species are included;
there is no Sim-photon-only selection. `SimPart_pdgId` stores DELPHI mass codes
in these samples, not PDG IDs. The link semantics come from the
[NanoAOD producer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L892).

The three curves use identical **gen energy** or **gen cos(theta)** bins:

```text
Gen → Sim    = N(S) / N(G)       saved-lineage coverage
Reco | Sim   = N(R) / N(S)       conditional truth-associated reconstruction
Gen → Reco   = N(R) / N(G)       total truth-associated reconstruction
```

A gen photon with several reco fragments counts once. The code checks
`R ⊆ S ⊆ G` and central-value closure
`εGR = εGS × εR|S` in each populated bin. A zero denominator is undefined and
not plotted. Each ratio has its own 68.27% Clopper–Pearson interval; the curves
are correlated, so their uncertainty bands must not be multiplied as independent
measurements. Using Sim or reco energy instead would introduce migration and
would not have this bin-by-bin closure.

**Gen → Sim is saved-lineage coverage, not detector transport or acceptance.**
A saved input Sim root can exist even when no detector response is recorded.
Conversely, missing stored lineage alone cannot identify transport thresholds,
geometry, collection pruning, or an association failure as the cause.
Tracing descendants can recover a reco origin, but cannot invent a missing Gen
anchor. The writer exports all SKELANA `NVECMC` entries without an additional
energy/angular selection; the completeness of the upstream simulation record
is not established by that fact.

### What the ROOT data show

An independent full-data audit finds that **all 1,795,749 stable gen photons
have exactly one direct saved Sim root**. `GenPart_simIdx` and `SimPart_genIdx`
are exact inverses, with zero missing links or disagreements. Root Sim
Px/Py/Pz/E are bit-identical to Gen in every case. These are saved creation/input
four-vectors, not calorimeter deposits or energy after detector interactions.
Here, a "Sim root" means the directly Gen-linked anchor of that photon's lineage,
not necessarily a root of the entire Sim graph: 138,318 Ztautau photon anchors
have Sim parents. The origin rule stops at that first Gen anchor.
Consequently, `Gen → Sim = 100%` and `Reco | Sim = Gen → Reco` in every
populated bin. The overlapping curves are intentional.

All 899,998 Beam ISR photons have a saved root, **zero saved descendants, and
zero associated reco Photon candidates**. This does not prove that they were
never transported: the saved Sim schema has no transport flag, detector hits,
or deposited-energy field. Reference simulation title files contain thresholds,
but sample-specific titles/logs are unavailable, so no loss is assigned to them.
`SimVtx_position` is not used to infer interaction locations; its producer
coordinate indexing is inconsistent with the DELPHI vertex layout. The integer
ancestry links used here pass the independent graph audit.

| Sample | Stable gen photons | Gen → Reco: all stable | Gen → Reco: all ISR | Gen → Reco: Non-beam ISR |
| --- | ---: | ---: | ---: | ---: |
| Zee | 371,987 | 3.389% | 0.585% | 2.175% |
| Zmumu | 313,532 | 3.713% | 0.557% | 2.073% |
| Ztautau | 505,036 | 16.202% | 0.581% | 2.160% |
| ZKK | 295,891 | 3.456% | 0.556% | 2.075% |
| Zpipi | 309,303 | 3.718% | 0.559% | 2.096% |

There are 6,983 unique reconstructed gen ISR photons, versus 6,757 with only
an immediate reco→Sim→Gen link. Sim ancestry recovers 226 additional gen ISR.
There are 7,020 ISR-associated reco candidates; 34 gen ISR photons have multiple
candidates, with a maximum of three. ISR with saved Sim descendants numbers
3,935 / 3,848 / 3,879 / 3,908 / 3,795 in sample order. This is a saved-lineage
diagnostic, not a calorimeter-hit or isolated-photon detection requirement.

The reco `Photon_fourMomentum` collection contains neutral electromagnetic
calorimeter candidates, not just gen photons. **362,929 of 573,796 reco candidates
have no stored Sim link.** They remain unassociated in nominal efficiency;
their physical origin is not established from truth. Thus Gen→Reco measures
reconstruction identifiable from the stored associations, not every possible
physical ISR contribution to calorimeter candidates.

### Reading the figures

Each group has these plots:

- `gen_lineage_*`: Gen / With Sim lineage / With reco, all in **gen coordinates**.
  Multiplicity counts the unique gen photons in each subset per event.
- `gen_sim_reco_*`: actual Gen / Saved Sim roots / Associated reco coordinates.
  The first group additionally shows **all reco**, including unassociated
  candidates, as black points. Associated reco may include several fragments
  from one gen photon. Gen and saved Sim curves coincide in these samples.
- `stage_efficiency_vs_energy` and `stage_efficiency_vs_cos_theta`: the three
  gen-unit efficiencies with Clopper–Pearson error bars.
- `lineage_multiplicity`: saved Sim objects (root plus descendants of all species)
  and reco candidates per gen photon, including zero reco. The y axis is a gen
  fraction. Sim ancestor and descendant energies are never summed.

Coordinate histograms use count/bin/event and sqrt(N) errors; event multiplicity
uses event fraction. There are no plot titles. CMS styling, default font sizes,
11×9 figures, DELPHI Simulation, and bold sample legends are retained.
Plots with log axes or color scales have fully linear `_linear.png` siblings.
Coordinate log-energy distributions have 60 bins; energy efficiency has roughly
three bins per decade to reduce sparse-bin fluctuations. Linear energy plots
use 1 GeV bins. All plotted ranges include every selected entry. Different bin
widths mean count/bin/event heights should be compared within the same scale.

The ISR group also retains the truth energy response, summed-reco energy residual,
opening angle, ancestry depth, matched-gen energy–angle maps, and generator
`x_gamma`/`|cos(theta)|` checks. Pair response counts every associated reco;
summed-reco response counts each matched gen once. Large positive energy tails
are retained; a single stored association does not establish their cause.
`sqrt(s)` for `x_gamma = 2 E_gamma / sqrt(s)` comes from the two Gen beams:
91.186996 GeV, versus metadata `Event_cmEnergy = 91.25` GeV.
Generator spectra show qualitative soft/forward enhancement consistent with
[QED shower emission](https://arxiv.org/abs/1410.3012); Beam ISR includes
[exactly collinear residual-radiation representatives](https://pythia.org/latest-manual/htmldoc/PDFSelection.html).
An independent same-process reference such as [KKMC](https://arxiv.org/abs/hep-ph/9912214)
is needed for precision agreement. The local Zee sample is s-channel annihilation.

### Auxiliary angular matching

`04_geometric_matching` retains direct truth, truth ancestry, angle-only, and
truth-plus-geometry comparisons. Angle-only uses closest-first one-to-one 3D
opening-angle matching of all stable gen photons to all reco candidates at
0.03 rad. Geometry fallback keeps truth matches and uses only reco with no
resolved Gen origin. These diagnostic matches never enter the nominal stage
efficiencies.

At 0.03 rad, angle-only ISR pairs have 6,661 correct origins, 18 different
origins, and 512 unknown origins. Known-pair purity is 99.73%, excluding unknown
origins; truth-pair acceptance is 6,691/7,020 = 95.31%. These have different
denominators. Truth plus geometry adds 479 diagnostic gen matches.

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
