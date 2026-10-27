# DELPHI Analysis

ROOT validation and BDT event classification for DELPHI NanoAOD.

## Setup

Use Python 3.12 and run commands from the repository root:

```bash
uv sync --locked
```

```text
scripts/
  01_sample_preparation/      download, conversion and train/val/test datasets
  02_sample_validation/      MC branch, track, gamma and RICH validation
  03_data_mc_comparison/     collision Data and MC comparisons
  04_event_classification/    BDT training, evaluation and inference
src/                         reusable Python code
config/pid/                  feature sets and hyperparameters
docs/                        data definitions and branch checks
data/                        source samples
output/{dataset,plots,models}/
```

Bash recipes contain defaults; append arguments to override them. Packages are
installed from `src/`, so no `PYTHONPATH` setting is needed. Data, outputs
and caches are excluded from Git.

Use the SQLite DB from the official [delphi-metadata v0.1.0 release](https://gitlab.cern.ch/liko/delphi-metadata/-/releases/v0.1.0),
placed at the repository root as `delphi-metadata.sqlite3` (excluded from Git).
The current local DB is already in place. Query it directly with the official CLI:

```bash
uv run dmd runquality breakdown --year 1994
uv run dmd runquality legend
```

These queries are independent of collision files. Once conversion finishes,
check the actual events and save `output/dataset/<campaign>/run_quality_runs.json`:

```bash
uv run python -B -m sample_preparation.run_quality
```

The checker joins native events to raw `fileSeq` by `(run, event)`. Its current
all-region reference requires TPC, VD, MUB and MUF in `[6,7]`, following the email.
TPC uses `tpc_0/1`, VD `mvx_a/c`, MUB `mub_b/d`, and MUF `muf_a/c`.
This reference requires both muon systems; the two-track comparison below
applies track-angle-dependent requirements.

Compare tracking-only, tracking+MUB, tracking+MUF and tracking+both using the
full 1994 metadata (independent of the collision conversion):

```bash
uv run python -B -m sample_preparation.run_quality --plot
```

The two plots in `output/plots/1994/run_quality/` show STIC luminosity retained
by each requirement and the MUB/MUF pass/fail overlap within tracking-selected
metadata. They measure available luminosity, not detector reconstruction efficiency
or the luminosity of the downloaded event subset. Every required flag must be
6 or 7; TPC/VD are not silently restricted to 7 alone.

## Prepare datasets

Pair native/raw events and write one train/val/test dataset per sample:

```bash
bash scripts/01_sample_preparation/prepare_dataset.sh \
    --input data/20260828_florian --output output/dataset/20260828_florian
```

Defaults are 60,000 training, 15,000 validation and up to 15,000 remaining test
events per sample. Explicit `--counts TRAIN VAL TEST` requires those exact sizes.
Sources are preserved and existing datasets are not overwritten.
See [data rules](docs/data.md) for the 14,999-event Florian Zmumu test split.

```text
output/dataset/<campaign>/
  metadata.json
  <sample>_<FSR>/{train,val,test}.root
    analysisLevel/Events
    raw_sdst/Events
```

The two trees contain the same events in the same order, keyed by
`(job_id, run, event)`. Preparation uses uproot and needs no ROOT installation.

## Validate samples

```bash
bash scripts/02_sample_validation/plot.sh
```

The default sample is ZKK PHOTOS-FSR. Select others with `--samples`, for example:

```bash
bash scripts/02_sample_validation/plot.sh --samples Zee Zmumu
```

Available samples are `ZKK`, `Zpipi`, `Zee`, `Zmumu` and `Ztautau`.
Plots use `train` by default; pass `--split val` or `--split test` to change it.
Each study's `plot.sh` runs all stages; numbered recipes run one stage.
Outputs follow the script names:

```text
output/plots/<campaign>/02_sample_validation/<split>/
  {00_overall,01_track,02_gamma,03_rich}/<numbered_stage>/
```

Both campaigns use a 1° gamma matching cone. See [photon definitions](docs/photons.md)
for selections and efficiency denominators. Old Florian RICH associations
support storage checks only; see [converter changes](docs/converter.md).
Gamma distributions are split into `01-01_sample_distribution_all` and
`01-02_sample_distribution_fiducial`. Fiducial applies the existing Gen EM
angular acceptance; the Reco distributions remain inclusive in both.
Gamma matching validation and results each have `all/` and `fiducial/` directories.
Both exclude Beam ISR; only `fiducial` applies the Gen EM angular acceptance.

## Compare Data and MC

```bash
bash scripts/03_data_mc_comparison/plot.sh
```

This study defaults to original Florian PHOTOS-FSR Zmumu and Zee MC. The MC-only
`02_sample_validation` stage still defaults to ZKK. `--samples` selects MC
templates; it does not label Data or choose its control channel. Both Data
selections run independently; any overlap is retained. Common plots include all
requested MC samples. Channel plots compare Data only with the corresponding
MC: Zmumu in `01_mumu`, Zee in `02_ee`. The other channel was a background
template diagnostic and is omitted from these channel comparisons.

| Stage | Requirement |
| --- | --- |
| Common quality | Data TPC/VD quality 6 or 7, before any track selection |
| Common tracks | Exactly two charged Part rows, `q1*q2 == -1`; finite positive momentum |
| Common backtoback | `pi - delta_phi < 0.35` radians |
| Muon channel | Angular acceptance → both standard muon IDs → regional muon quality |
| Electron channel | Angular acceptance → both standard electron IDs → regional EM quality |
| Final mass window | `70 < m_ll < 110` GeV, after regional quality in each channel |

| Channel / region | Theta (degrees) | Required detector quality |
| --- | --- | --- |
| Muon barrel | 53–88.5 or 91.5–127, inclusive | `mub_b` and `mub_d` |
| Muon endcap C / A | 20–42 / 138–160, inclusive | `muf_c` / `muf_a` |
| Electron barrel | 44–88 or 92–136, inclusive | `hpc_0` and `hpc_1` |
| Electron forward C / A | 20–35 / 145–160, inclusive | `emf_c` / `emf_a` |

TPC/VD always uses `tpc_0/1` and `mvx_a/c`. A is negative z; C is positive z.
Every track contributes its needed flags, including mixed barrel/forward pairs.
Electron selection has no MUB/MUF condition; muon selection has no HPC/FEMC condition.
Both channels apply regional quality after ID and before the final mass
selection. Quality 6 or 7 is required only for Data and only for the regions
crossed by the two tracks. Muon requirements are:

| Muon pair | Required quality |
| --- | --- |
| Barrel–Barrel | `mub_b` and `mub_d` |
| Endcap–Endcap | The visited `muf_a` / `muf_c` sides |
| Barrel–Endcap | `mub_b`, `mub_d` and the visited MUF side |

Mixed pairs remain in the quality logic and CSV counts. Mass plots combine
the accepted regions.

HPC 0/1 are both required because their numbering has not been established as an
A/C split. Electron angles are a conservative **working choice**, not a selection
specified in the supplied slides. The forward intervals stay inside both tracking
coverage and the published FEMC geometry. Muon angles follow
[CERN-PPE-95-194, section 16](https://repository.cern/records/kz87a-w1k25/files/ppe-95-194.pdf?download=1).

Native PID rows are joined by `Muid_partIdx` / `Elid_partIdx`. Standard Muon ID is
`(Muid_tag & 4) != 0`; `Muid_tag > 0` would also accept HCAL-only tags. Standard
Electron ID is `Elid_tag >= 4`. Missing or ambiguous PID links stay unknown and
fail the ID stage. Raw PID fields are not used. Channel mass distributions and
final mass cuts both use track momenta with the corresponding muon or electron
mass hypothesis. The stored-four-vector mass in common plots is diagnostic only.
The standard muon bit uses chamber hits/refits; the HCAL-only bit is separate.
Barrel ElID can also use RICH/OD and bremsstrahlung information through multiple
OR paths, which the final native tag does not identify. TPC+EM quality does not
establish that every upstream ElID input is fully qualified.
The electron ID stage is therefore a stored-tag enrichment/response check;
quantitative ID validation needs additional quality and independent-tag studies.

The two-track control study applies its mass window last within each channel.
The azimuthal 0.35-radian cut is **not** Sanghyun's 3D acollinearity selection.
`Back-to-back (phi)` means `pi - delta_phi < 0.35` rad: the tracks are within
about 20 degrees of opposite azimuthal directions. `Both muon IDs` / `Both
electron IDs` means that both stored native tags pass the standard-ID definitions
above. This is an extra control-sample selection after common cuts and
angular acceptance, not a BDT cut. Both channels use acceptance → standard ID
→ regional quality, followed by the final mass selection.
Sanghyun's momentum-resolution and tag-and-probe references (slides p.7 and p.10) use
different momentum, angle and quality requirements, which are not silently mixed
into this baseline. Dietrich's metadata slides p.10–12 motivate detector-specific
quality; our agreed baseline remains 6–7, rather than a universal strict-7 cut.

There is no BDT selection, photon veto, `nPart == 2`, 35-GeV momentum/energy cut,
or electron veto in this baseline. Common counts include zero-ID candidates. The
standard-ID sample is enriched and selection-biased; it is not an unbiased PID
efficiency measurement. Independent tag-and-probe and background subtraction
remain separate studies.

The local Zee generation configuration enables `WeakSingleBoson:ffbar2gmZ`
with electron decays, consistent with all 90,000 saved Zee Gen records. Original
per-job generation logs are unavailable; full Bhabha coverage is not established.
Treat this MC as an s-channel gamma*/Z template. Forward
dielectron Data can therefore differ through missing t-channel/interference
physics as well as detector effects. A full Bhabha sample is needed before
interpreting those shape differences quantitatively. See the
[Pythia process definitions](https://gitlab.com/Pythia8/releases/-/raw/pythia8313/share/Pythia8/xmldoc/ElectroweakProcesses.xml).

Data defaults to `data/20260828_florian/20261026_short94_c2`; MC defaults to
`data/20260828_florian`. Only completed native/raw conversions with exact 1:1
physics keys enter. Run Quality is applied to Data only, keyed by `(run, fileSeq)`.
TPC/VD quality is evaluated on all native events first, preserving original row
indices for event and PID links. Exactly-two is then counted before the OS stage
rejects pairs with invalid charge or non-finite/zero momentum; no third track is
dropped to promote an event.

Outputs: `output/plots/<MC campaign>/03_data_mc_comparison/two_track/`:

```text
00_common/
  00_cutflow.png
  01_all_charged_parts_{p,cos_theta,inv_mass}.png
  02_tracking_quality_*.png
  03_two_charged_tracks_*.png
  04_opposite_sign_*.png
  05_backtoback_*.png
  06_selected_*.png
01_mumu/
  00_cutflow.png
  01_common_*.png
  02_angular_acceptance_*.png
  03_muon_id_*.png
  04_detector_quality_*.png
  05_mass_window_*.png
  06_selected_*.png
02_ee/
  00_cutflow.png
  01_common_*.png
  02_angular_acceptance_*.png
  03_electron_id_*.png
  04_detector_quality_*.png
  05_mass_window_*.png
  06_selected_*.png
cutflow.csv                    common/channel/region event counts; MC quality is N/A
inputs.csv                     input snapshot and invalid/unknown counts
histogram_flows.csv             entries outside display ranges or non-finite
```

Every file has one axis and shows the sample **after the named cumulative
selection**. Read the numbered files in order. The upper left shows only the new
cut at that stage; the channel baseline is labelled `Common`. Common plots
include the requested MC samples; channel plots use Data and the matching
MC only. Cutflows use Cut on x and Events on y. Channel cutflows start with the
full No cut count, followed by Common and the channel's selection order above. At
Data-only quality steps, MC counts stay constant. Plot labels omit the `(Data)`
suffix; quality selections still apply only to Data. The common cutflow uses
log-y.

`06_selected_*.png` repeats each directory's final distributions without the cut
annotation, keeping the same axes and normalization: backtoback in `00_common`,
and the mass-window selection in `01_mumu` and `02_ee`.

`00_common` progresses through No cut → TPC/VD quality → exactly two charged
Parts → OS (including pair validity) → backtoback.
The backtoback cut remains, but its separate variable-distribution plots are not
produced. Momentum and cos(theta) use all charged Parts remaining at each stage.
The invariant mass
uses the **two highest-momentum charged Parts** before the exactly-two cut,
summing their stored `Part_fourMomentum` 4-vectors. Once exactly two charged
Parts are required, it uses that pair. This mass is diagnostic and does not
select events. With fewer than two charged Parts, mass is undefined, not zero. These
entries, negative mass-squared and other undefined observables are recorded as
non-finite in `histogram_flows.csv`, without an extra plotting cut.

Both channels start from Common without a mass cut, then use angular acceptance
→ both stored standard IDs → regional quality. Electron quality requires
`hpc_0` and `hpc_1` for barrel tracks and the visited `emf_a` / `emf_c` sides
for forward tracks. Each channel then applies `70 < m_ll < 110` GeV using the
same fixed-mass `m_ll` shown in its plots. Channel mass plots retain the full
0–140 GeV display.

Cutflows show raw event counts; distribution curves are independently normalized
to compare shapes. Electron barrel `E_HPC/p` uses `Part_hpcTotalShowerEnergy` on
log-y, including zero values. No charged FEMC energy branch is available, so
`Photon_emEnergy` is not substituted. No dE/dx/RICH measurement-validity cuts
are introduced.

The 0–140 GeV mass and 0–70 GeV momentum ranges are display ranges, not cuts.
Over/underflows are recorded and retained in normalization; non-finite observable
entries are counted separately. Each MC curve is independently normalized,
including overflow, and is a shape template rather than a weighted prediction.

## Train and evaluate

Choose a campaign and one of `detector`, `detector_pid` (default),
`detector_rich` or `detector_pid_rich`. All use the `default` hyperparameters:

```bash
export CAMPAIGN=20261015_jongwon FEATURE_SET=detector_pid
bash scripts/04_event_classification/01_model_train.sh
bash scripts/04_event_classification/02_model_evaluate.sh
bash scripts/04_event_classification/03_model_apply.sh
```

Training fits `train` and uses `val` for early stopping. Use `test` for the final comparison.
See [BDT inputs](docs/bdt.md) for features, classes and missing-value rules.

| Output | Location |
| --- | --- |
| Models and metrics | `output/models/<campaign>_<feature>_default/` |
| Evaluation plots | `output/plots/<campaign>/04_event_classification/02_model_evaluate/<feature>_default/` |
| Predictions | `output/dataset/<campaign>/predictions/<feature>_default/` |

Experiment summaries are stored in `output/models/results.csv`.
See the [documentation index](docs/README.md) for branch checks and known limits.
