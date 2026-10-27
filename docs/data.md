# Data

## Files and campaigns

| Source file | Dataset tree | Contents |
| --- | --- | --- |
| `nanoaod.root:Events` | `analysisLevel/Events` | Reconstructed Part/Photon objects and Gen–Sim–Reco links |
| `nanoaod_raw_sdst.root:Events` | `raw_sdst/Events` | Detector measurements, raw tracks/PID and calorimeter clusters |

`20260828_florian` uses the old converter. `20261015_jongwon` includes the
[Gamma/RICH fixes](converter.md). Existing files need reconversion to receive fixes.

Preparation pairs events by `(job_id, run, event)`, removes raw setup records,
excludes all rows with duplicate event keys, preserves branch structures, and adds `job_id`.
The paired trees have identical event order; train/val/test do not overlap.
Source files are retained. `metadata.json` records input paths, excluded jobs/rows
and actual split sizes.

| Campaign / sample | Train | Val | Test |
| --- | ---: | ---: | ---: |
| Jongwon, all five samples | 60,000 | 15,000 | 15,000 |
| Florian, except Zmumu | 60,000 | 15,000 | 15,000 |
| Florian Zmumu | 60,000 | 15,000 | 14,999 |

Florian Zmumu is missing one original event: job 16, run −3016, event 1271.
The default preparation handles this without a sample-specific condition:

```bash
bash scripts/01_sample_preparation/prepare_dataset.sh
```

It assigns 60,000 training and 15,000 validation events, then up to 15,000
remaining test events. Explicit `--counts TRAIN VAL TEST` still requires exact sizes.
Preparation defaults to ZKK PHOTOS-FSR. Use `--samples Zee Zmumu` to select other
PHOTOS-FSR channels, or list all five to prepare them together.
Dataset directories retain the FSR suffix.

## Collision data

The converted 94C2 files are under
`data/20260828_florian/20261026_short94_c2/<file>/`, alongside the Florian MC
campaign. Each directory contains the original `.al`, `nanoaod.root`,
`nanoaod_raw_sdst.root` and `conversion.log`.

Download, conversion and event Run Quality commands default to this location.
Data/MC plots use `scripts/03_data_mc_comparison/plot.sh`; MC-only branch and
matching validation remains in `scripts/02_sample_validation/`.

The conversion recipe uses the prepared local Docker image `delphi-analysis:local`,
`.cache/short94_c2/{src,build,validate.py,provenance.txt}` and
`.cache/sample/condition-data`. These assets are local prerequisites; the recipe
does not build or download them.

## Join objects by index

- Original source rows must be paired by event keys, not file row number.
- PA IDs and native Part/VECP rows are different index spaces. Corrected
  `Part_paIdx` and `HaidRaw_paIdx` are zero-based, with `−1` for no link. PA 0 is valid.
- Old Florian `HaidRaw_paIdx` contains one-based VECP numbers. Subtracting one
  does not turn it into a raw-track PA ID.
- Native Trac/HAID/RICH/PID arrays follow Part order and include neutral placeholders.
  Their lengths are not identified-particle counts.
- Raw collections have their own rows; use stored PA/track indices to associate them.
- Native `Vtx_incomingIdx/outgoingIdx` are one-based with 0 for no link.
  `SimVtx_*` indices are zero-based with `−1` for no link.
- Native track refits and B-tag primary-vertex replacement can make track/Vtx
  values differ from raw banks. See [classifier input checks](branches/classification.md).

## Interpret stored values

| Field | Meaning |
| --- | --- |
| `GenPart_pdgId` | Standard PDG code |
| Native `Part_pdgId`, `SimPart_pdgId` | DELPHI mass code: gamma 21, electron 2, positron −2 |
| `Photon_emEnergy` | FEMC energy; HPC is stored separately |
| `Stic_directionFromMain` | Momentum in GeV, not a unit direction |
| `HadShower_direction` | XYZ impact position in cm, not a unit direction |
| `Event_totalChargedEnergy` | Charged energy after Team4 selection |
| `Event_totalEMEnergy` / `Event_totalHadronicEnergy` | Neutral EM / neutral hadronic energy |

Zero, `−1` and empty arrays have field-specific meanings. Finite values can
still come from incorrect words. Negative/extreme track χ² and momentum
outliers need upstream investigation; a correct copy is not a fit-quality check.
Raw MUID fields and native composite perigee/matrix fields still have known
porting defects; see the [native](branches/native.md), [raw](branches/raw.md)
and [current classifier](branches/classification.md) checks.

Use Gen ancestry for ISR classification. Time/lifetime information is absent
in these DST samples. Use unique, consistent stored Gen origins; leave
conflicting or missing links unresolved and stop Sim ancestry at the first
valid Gen anchor. Photon candidates are EM candidates, not confirmed PID photons;
see [photon definitions](photons.md).

Florian RICH has invalid layout/association information. Jongwon fixes storage
and PA links; RICH quality cuts, calibration and PID performance are separate checks.
In the [scan tables](README.md), `values` counts flattened numbers and
`subnormal` counts nonzero floats below the smallest normal magnitude.
