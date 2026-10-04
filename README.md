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
`plot_isr_photons.py` makes the comparison and efficiency figures,
`plot_2d.py` makes matched-ISR maps, and `plot_response.py` and
`plot_generator.py` contain response and generator checks.

Stable gen photons have `GenPart_pdgId == 22` and `GenPart_status == 1`;
reco photon candidates use `Photon_fourMomentum`. The `Photon` collection
contains neutral electromagnetic calorimeter candidates, not just gen-linked
photons. Gen ancestry divides stable photons into Non-beam ISR, Beam ISR, FSR,
and Others. Beam ISR has zero transverse momentum and a direct parent from the
incoming parentless electron or positron. Non-beam ISR is the remaining ISR.
The exact ancestry rules are in `data.py`. There is no gen energy threshold,
fiducial-angle cut, or matching energy cut. Invalid photon four-vectors stop
the run rather than being removed.

Stage 01 stacks the four stable-gen origin categories in the 1D energy,
cos(theta), and phi distributions. Stage 02 removes Beam ISR from the gen
stack and gen total. Both stages overlay the same set of **all** reco photon
candidates as black points with Poisson error bars. Corresponding figures in
both stages use the same bins, with gen and reco divided by the same number of
events. This compares populations; reco candidates need not match the gen
photons underneath them. Multiplicity shows the separate gen-origin event
distributions, the gen total, and reco points. The origin distributions cannot
be stacked into the total multiplicity distribution.
Photon-coordinate histograms show `Nγ / bin per event` (bin count divided by
the number of events), while multiplicity uses `Event fraction`. A photon with
undefined phi is displayed in the phi=0 bin solely for plotting.

The nominal result is **truth-associated ISR reconstruction efficiency**:
the fraction of gen ISR photons with at least one associated reco candidate.
Follow `Photon_partIdx → Part_simIdx`. If `SimPart_genIdx` is negative, follow
the Sim parent through `SimPart_originVtxIdx → SimVtx_incomingIdx`, repeating
until the first valid Gen index. That first Gen particle must be a stable ISR
photon. Direct links have ancestry depth zero; conversion and shower
descendants can have positive depth. The link semantics are defined by the
[NanoAOD producer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L892).
`SimPart_pdgId` stores DELPHI mass codes in these files, so the study does not
apply a PDG selection to it.

Nominal matching has **no opening-angle or energy requirement**. A gen photon
with multiple linked reco candidates counts once in the efficiency numerator.
Reco candidates with no Sim association remain unassociated. The result thus
measures reconstruction identifiable from the stored truth graph; it does not
prove that every unassociated candidate is unrelated to ISR.
Stage 04 uses all ISR and Non-beam ISR as separate denominators, with 68.27%
Clopper–Pearson intervals. Both include detector acceptance.

Stage 03 plots nominal matched-gen photons and truth response. Energy response,
opening angle, and ancestry depth count every associated reco candidate;
`isr_reco_descendants` counts each gen ISR once, including unmatched zeroes.
The summed-energy residual uses the sum of all associated reco energies for
each matched gen ISR. These candidate energies need not obey isolated-photon
energy conservation because the collection and its single stored truth link
can represent shower fragments and merged calorimeter deposits. Response
plots retain their complete observed ranges.

Angular matching is an auxiliary study in `05_geometric_matching/`:

- **Direct truth:** only the immediate `SimPart_genIdx` link.
- **Truth ancestry:** the nominal Sim-parent tracing above.
- **Angle only:** closest-first, one-to-one matching of all stable gen photons
  to all reco candidates, ignoring truth, with 3D opening angle < 0.03 rad.
- **Truth + geometry:** keep nominal truth matches, then match unmatched gen
  photons to reco candidates with no resolved Gen origin using the same angle
  rule. Every candidate with a known Gen origin is excluded from this fallback.

The angle-cut validation shows **known-pair purity**, correct-origin pairs
divided by correct plus different-origin pairs among angle-selected ISR pairs
with known truth. It also shows **truth-pair acceptance**, the fraction of all
truth-associated reco–ISR pairs inside the cut. Unassociated pairs are shown
separately in the count scan and are not assumed correct. These two fractions
have different denominators. The nominal efficiency never uses the angle cut.

Generator checks in stage 01 show `x_gamma = 2 E_gamma / sqrt(s)` and
`|cos(theta)|`, split into Non-beam and Beam ISR. Their vertical axes are
bin-width-normalized densities per event. `sqrt(s)` comes from the first two
generator beam four-vectors, 91.186996 GeV in these samples; `Event_cmEnergy`
instead stores 91.25 GeV and is retained only as metadata in the summary.
The generator checks test the qualitative soft and forward enhancement of
[QED shower emission](https://arxiv.org/abs/1410.3012). Beam ISR contains the
[exactly collinear residual-radiation representatives](https://pythia.org/latest-manual/htmldoc/PDFSelection.html),
so its photon count is not a resolved-emission prediction. No arbitrary `1/E`
theory curve is overlaid. Precision agreement would require an independent
reference such as [KKMC](https://arxiv.org/abs/hep-ph/9912214) with the same
process and selections. The local Zee configuration is s-channel annihilation,
not a full Bhabha sample.

Stages 01 and 02 use the same scales: `gen_reco_energy.png` has log energy and
count axes, while the cos(theta) and phi defaults have log count axes. Their
`_linear.png` counterparts have only linear axes; for example,
`02_gen_reco_wo_beamISR/Zee/gen_reco_energy_linear.png`. Multiplicity is already
linear in both stages. Other figures with a log axis or color scale also have a
fully linear `_linear.png` version. Linear energy plots use 1 GeV bins for 1D
distributions and efficiency; matched energy–angle maps use 2 GeV bins, and
energy-response maps use 1 GeV bins. Their ranges include every plotted
photon. Linear efficiency plots show 0–1.05 so the
Clopper–Pearson intervals in sparse high-energy bins remain visible. Since
energy bin widths differ, `Nγ / bin per event` heights should be compared only
within the same scale.

Plots are grouped under `plots/20260828_florian/isr_photons/`. Each sample has
one axis per 1D figure. Matched-gen 2D energy–cos(theta) maps in stage 03 use
common bins and color scales across samples.

```text
01_gen_reco/<sample>/              gen/reco 1D overlays
02_gen_reco_wo_beamISR/<sample>/  gen/reco 1D overlays without Beam ISR in gen
03_matching/<sample>/            nominal matched gen ISR, Sim ancestry, and response
04_efficiency/<sample>/          nominal all-ISR and Non-beam-ISR efficiencies
05_geometric_matching/<sample>/  auxiliary matching comparisons
05_geometric_matching/           pooled angle-cut counts and truth validation
study_summary.json              counts, integrated efficiencies, and validation metrics
```

The full run contains 449,999 events. The nominal counts agree with an
independent traversal of the ROOT Sim graph:

| Sample | Direct matched ISR | Ancestry matched ISR | All-ISR efficiency | Non-beam-ISR efficiency |
| --- | ---: | ---: | ---: | ---: |
| Zee | 1,395 | 1,441 | 0.585% | 2.175% |
| Zmumu | 1,320 | 1,372 | 0.557% | 2.073% |
| Ztautau | 1,377 | 1,430 | 0.581% | 2.160% |
| ZKK | 1,329 | 1,367 | 0.556% | 2.075% |
| Zpipi | 1,336 | 1,373 | 0.559% | 2.096% |

Truth ancestry adds 226 unique gen ISR photons (6,757 → 6,983). Beam ISR has
zero truth-associated reconstruction in every sample. There are 7,020 linked
reco candidates: 34 matched gen ISR photons have multiple candidates, with a
maximum of three. For example, ZKK `job_0`, entry 1637 has two reco candidates
whose Sim-parent paths both reach gen ISR index 6; this contributes one to the
nominal numerator.

At 0.03 rad the angle-only study has 6,661 correct-origin, 18 different-origin,
and 512 unassociated ISR pairs. Its known-pair purity is 99.73%, excluding the
512 unassociated pairs. Truth-pair acceptance is 6,691/7,020 = 95.31%.
Truth plus geometry adds 479 diagnostic gen matches; none enters the nominal
efficiency. The stored Sim link is absent for 362,929 of 573,796 reco
candidates, so missing truth remains a material limitation.

In Zmumu, 92.89% of Non-beam ISR photons have energy below 1 GeV and 54.15%
have `|cos(theta)| > 0.9`, consistent with qualitative soft/forward enhancement.
These are reported fractions, not selection cuts. The median pair energy
residual is about −3% across the samples, with large positive tails retained
in the response figures. A single stored association does not establish the
cause of those tails.

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
