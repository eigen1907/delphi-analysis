# BDT

Use the same feature set and hyperparameters to compare Florian and Jongwon.
Train on `train`, use `val` for early stopping, and compare final results on `test`.

Class IDs are `Zee=0`, `Zmumu=1`, `ZKK=2`, `Zpipi=3`, `Ztautau=4`.
Confusion matrices display `Zee`, `Zmumu`, `Ztautau`, `ZKK`, `Zpipi`; this does not change the IDs.

## Inputs

| Feature set | RAW inputs | Model columns |
| --- | --- | ---: |
| [detector](../config/pid/features/detector.json) | Detector collections below | 743 |
| [detector_pid](../config/pid/features/detector_pid.json) | Detector + `MuidRaw`, `ElidRaw`, `HaidRaw` PID | 856 |
| [detector_rich](../config/pid/features/detector_rich.json) | Detector + eight RICH measurements | 823 |
| [detector_pid_rich](../config/pid/features/detector_pid_rich.json) | Detector + PID + eight RICH measurements | 936 |

The default is `detector_pid`. Counts include the five analysis-level inputs below.
The JSON files define the full branch lists; the saved model's `features` records the column order.

Detector collections from `raw_sdst/Events`:
`Event`, `TracRaw`, `EmShower`, `EmLayer`, `HadShower`, `HadHit`,
`Stic`, `MtpcRaw`, `VdAssocHit`, `VdUnassocHit`, `Vtx`.

The eight RICH measurements are:

```text
HaidRaw_thetaGas  HaidRaw_sigmaGas  HaidRaw_nphGas  HaidRaw_nepGas
HaidRaw_thetaLiq  HaidRaw_sigmaLiq  HaidRaw_nphLiq  HaidRaw_nepLiq
```

`detector_rich` adds only these eight branches. `detector_pid_rich` also retains
the PID tags and `nHaidRaw`. Neither adds RICH quality or sector flags.

Every feature set adds these branches from `analysisLevel/Events`:

```text
Event_chargedMult       Event_neutralMult
Event_totalChargedEnergy
Event_totalEMEnergy     Event_totalHadronicEnergy
```

Gen/Sim truth, event/job IDs, and labels are excluded from model inputs.

## Extraction

- Keep ten tracks ordered by increasing `abs(TracRaw_invR)`. Keep the 15 covariance elements for the first four tracks.
- Join track-related collections by their stored PA or track indices. Do not join particle rows across trees by position.
- Pad missing slots with NaN. For old Florian data, unresolvable HaID associations remain NaN while configured event-level counts are retained. Extra RICH columns therefore provide no new measurements for Florian.
- Check that the two trees have identical event counts, job IDs, run numbers, and event numbers.

The [porting review](branches/classification.md) identifies errors in Raw MUID inputs.
Identical feature names alone do not guarantee identical information across campaigns.

## Training defaults

The only hyperparameter preset is [default.json](../config/pid/hyperparameters/default.json).

| Parameter | Value |
| --- | --- |
| `objective`, `eval_metric` | `multi:softprob`, `mlogloss` |
| `tree_method`, `grow_policy` | `hist`, `lossguide` |
| `n_estimators`, `early_stopping_rounds` | 70, 8 |
| `max_leaves`, `max_depth`, `max_bin` | 5, 0, 32 |
| `learning_rate`, `min_child_weight` | 0.25, 30 |
| `subsample`, `colsample_bytree` | 0.7, 0.4 |
| `reg_alpha`, `reg_lambda` | 0, 1 |
| `importance_type` | `gain` |
| Runtime seed, CPU threads | 1907, 32 |

`max_depth=0` means no depth limit; `max_leaves` limits tree size.
The training command overrides the preset's `n_jobs=-1` with `--ncpu` (default 32).

See the [README](../README.md) for recipes and the
[extractor](../src/event_classification/bdt.py),
[loader](../src/event_classification/dataset.py), and
[trainer](../src/event_classification/train.py) for implementation.
