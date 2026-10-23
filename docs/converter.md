# Converter fixes

Commits `ff6051a` (Gamma) and `47bc812` (RICH), after baseline `706d1dc`, change
the ROOT converter only. DELSIM, DELANA and the original events are unchanged.
The fixes apply to newly converted output; **existing `20260828_florian` files
retain the old values**.

Paired datasets rename the original `Events` trees to `analysisLevel/Events`
and `raw_sdst/Events`.

## Native changes

| Branch | Change |
| --- | --- |
| `Part_paIdx` | Added actual 0-based PA ID for each Part/VECP row, or `-1` if unavailable |
| `Rich_theg/sigg/nphg/nepg/flagg`, `Rich_thel/sigl/nphl/nepl/flagl` | Read the current event's decoded gas/liquid banks with correct COMMON and auxiliary types; preserve fractional expected photon counts |
| `Part_sticShowerEnergy`, `Part_sticShowerTheta`, `Part_sticShowerPhi`, `Part_sticNumTowers`, `Part_sticChargedTag`, `Part_sticSiliconVertexPos` | Correct STIC stride from 6 to 9 and silicon position from `KSTIC(6)` to `KSTIC(9)` |
| `Photon_sticShowerEnergy`, `Photon_sticNumTowers`, `Photon_sticChargedTag` | Apply the same STIC stride correction |
| `SimPart_genIdx`, `GenPart_simIdx` | Use Gen/LU ordering `ISTLU/ILUST`, replacing SH ordering `ISTSH/ISHST` |
| `SimVtx_masscode`, `SimVtx_position` | Read `KVTX(5)` and `QVTX(6..8)`, replacing slot 4 and slots 1..3 |
| `PhotonConv_simPhotonIdx` | Replace `KPHOC(3)-1` with a verified direct Sim-gamma link through the saved conversion vertex and reciprocal parent Part; otherwise store `-1` |

The STIC correction can change Photon membership and detector flags even though
the selection expression is unchanged. `Part_paIdx` is the only added native
branch; no existing branch was removed or renamed.

## Raw SDST changes

| Branch | Change |
| --- | --- |
| `HaidRaw_thetaGas/sigmaGas/nphGas/nepGas/flagGas`, `HaidRaw_thetaLiq/sigmaLiq/nphLiq/nepLiq/flagLiq` | Use the same decoded RICH banks as native; remove invalid/stale values and preserve fractional expected counts |
| `HaidRaw_paIdx` | Replace the 1-based VECP row with the actual 0-based PA ID, or `-1` |
| `HaidRaw_richQuality` | Widen int8 to int32 to preserve the full status word |
| `ElidRaw_tag`, `ElidRaw_gammaConvTag` | Read PA identification bits 7..9 and 10..12, replacing ELID header words |
| `ElidRaw_refitMomentum` | Read decoded `QELID(3..5)` |
| `EmShower_nLayers`, `EmShower_layerPattern` | Store `-1` for unavailable FEMC fields instead of reading beyond the record; retain HPC values |

The seven existing HAID PID tags (`kaonDedx`, `protonDedx`, `kaonRich`,
`protonRich`, `pionRich`, `kaonCombined`, `protonCombined`) were not recalculated.
The PA-index fix enables correct associations. FEMC was already skipped by the
`EmShower_nLayerEnergies`/`EmLayer` loop, so those values did not change.

## Definitions clarified without changing values

| Field | Meaning |
| --- | --- |
| Native `Part_pdgId`, `SimPart_pdgId` | DELPHI mass codes, not PDG: gamma=21, electron=2, positron=−2 |
| Native `Photon_emEnergy` | FEMC energy; HPC energy is separate |
| Native `Photon_sticChargedTag` | Current codes 0/1/2/3: no information / hard photon / soft photon / electron; legacy negative codes may occur |
| Native `PhotonConv_pairType` | Current positive codes 21–24; legacy codes may be negative |
| Native/raw Gen vertex | Position in cm; no usable time or lifetime information |
| Raw `Stic_directionFromMain` | MAIN momentum in GeV, not a unit direction |
| Raw `nHaidRaw` | VECP row count, including neutral placeholders |

Native `Haid_richQuality` was already int32. ROOT 6.36 API compatibility changes
do not alter branch names or meanings.

## Validation and limits

The original ZKK PHOTOS-FSR SDST and full FADANA samples were reconverted for
1,000 events each. A final 100-event conversion of each format matched every
column of the preceding corrected output.

| Check | Result |
| --- | --- |
| Native/raw RICH matched by PA | All 11 fields agree for 5,388 SDST rows and 5,090 full-FADANA rows; no subnormal values, duplicate PAs or missing track associations |
| Gamma sample | 1,280 Photon candidates and 68 conversions; 25 direct Sim-gamma links and 43 unresolved conversion origins |
| Sim vertices | 1,643 finite positions; 5,282 Gen–Sim position pairs agree |
| Raw ELID | Tags and momentum agree with 1,908 same-PA native rows; the other 177 raw columns are unchanged |
| Shower layers | 316 FEMC metadata entries are `-1`; 2,469 HPC showers and 12,775 layers retain their values |

- PA 0 is valid. Join PA and VECP/Part through explicit indices, not equal row numbers.
- RICH uses gas COMMON length 5 rather than 6, INTEGER NEP auxiliaries and REAL flag auxiliaries. `PSBEG` initializes the decoded banks; no extra reset was needed.
- RICH coverage means the fraction of tracks with a positive normal angle. It does not validate quality-bit cuts or PID performance.
- Photon is an EM-candidate collection. The converter does not invent missing STIC truth links or unresolved conversion origins; use saved parent/daughter ancestry.
- Empty raw EM/ELID in full FADANA does not mean no detector response: the raw EM reader expects SDST format.
- Electron refit momentum has a tail near 714 GeV that needs investigation before physics use.
- Native track offsets, raw MUID and packed TPC fields remain outside these commits. See the [native](branches/native.md), [raw](branches/raw.md) and [classification](branches/classification.md) reviews.

Validation files are retained in `data/20261014_stic_validation/`. See the
[project README](../README.md) for commands and [photon validation](photons.md)
for selection definitions.
