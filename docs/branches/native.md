# Native branches

This audit covers the original `20260828_florian` PHOTOS-FSR files,
`nanoaod.root:Events`. The scan on 2026-10-13 checked 100 files and 449,999 events:
90,000 per sample, except Zmumu with 89,999. It examined 310 numeric payloads
among 346 fields; the [scan CSV](../scans/native.csv) preserves the results.
Finite values can still come from the wrong memory word.

Paired datasets store this tree as `analysisLevel/Events`. The
[converter fixes](../converter.md) apply to newly converted files, **not to the
existing Florian files**. For model inputs, see the
[classification review](classification.md).

## Branch status

| Branch / collection | Florian finding | Current status / use |
| --- | --- | --- |
| `GenPart_*` | Momentum, PDG, status and ancestry checked; no time or lifetime information | Positions are in cm. Gen–Sim links were corrected by the converter fixes |
| `Part_*`, `Trac_partIdx` | Stored momentum and track-to-Part links are usable | Read native track momentum through `Trac_partIdx` and `Part_fourMomentum` |
| `Part_stic*`, `Photon_stic*` | COMMON stride was 6 instead of 9; silicon position read slot 6 | Fixed: stride 9 and silicon slot 9 |
| `Photon_*` | Neutral EM candidates; the STIC bug also affected candidate membership and detector flags | Fixed STIC values do not make this a photon-PID-selected collection |
| `PhotonConv_*` | All 45,543 direct Sim indices were `-1`; some daughter links were also missing | Use valid daughter/parent ancestry. The fixed converter writes direct Sim links only when verified |
| `Rich_theg/sigg/nphg/nepg/flagg` | Wrong gas COMMON length; all 1,699,997 nonzero `sigg` values were subnormal | Fixed converter reads the decoded gas bank |
| `Rich_thel/sigl/nphl/nepl/flagl` | All 2,099,957 nonzero `nepl` values were subnormal; flag type was wrong; other liquid values could be stale | Fixed converter reads the decoded liquid bank |
| `Haid_*`, `Haidc_*`, `Haide_*`, `Haidn_*`, `Haidr_*` | Separate PID banks; the Rich-array bugs do not establish that every PID tag is wrong | `Haid_richQuality` was already int32. `Haid_sign` is a PID-method word, not charge |
| `Trac_perigee`, `Trac_weightMatrix` | Read `QTRAC(3..7)` and `(8..22)` instead of `(4..8)` and `(9..23)` | Still outside the Gamma/RICH fixes. Do not reconstruct momentum or errors from these fields |
| `Trac_impParRPhi`, `Trac_impParZ` | Correct scalar slots 4 and 5 | Unaffected by the composite perigee error |
| `Trac_chi2NoVD`, `Trac_chi2VDHits` | 91,446 negative NoVD chi-squared values with positive ndf; large VD-hit tails | Also seen in raw data. Neither a normal chi-squared interpretation nor a universal failed-fit sentinel is established |
| `SimVtx_position`, `SimVtx_masscode` | Position contained index bits; mass code read the wrong slot | Fixed to `QVTX(6..8)` and `KVTX(5)`. Integer ancestry is a separate field |
| `Event_totalChargedEnergy`, `Event_totalEMEnergy` | Selected charged energy and neutral EM energy | They need not equal sums over all Part objects |
| `Elid`, `Muid`, `Btag`, `Uter`, `V0`, `Vtx` | This audit checked storage, not complete PID or fit performance | Check each branch definition and validity before use |

## Missing information

These are observations in Florian. An all-zero field alone does not prove a
converter defect.

| Branch | Limitation |
| --- | --- |
| `BeamSpot_errorFlag`, `BeamSpot_x/y/z`, `BeamSpot_sigmaX/Y/Z` | Remain at their initial zero values; error flag 0 does not establish a valid beam spot |
| `V0Hyp_kind/probability/vertex/momentum`, `V0Hyp_impactXY/impactXYErr/impactZ/impactZErr/impact3D/impact3DErr` | `nV0Hyp=8*nV0`, but all physical values are zero; SDST version 105 does not fill the old hypotheses |
| `Trac_resRPhiFirstPoint`, `Trac_errorResRPhiFirstPoint`, `Trac_resZFirstPoint`, `Trac_errorResZFirstPoint` | No residual information |
| `Dedx_width`, `Dedx_gapWires` | Explicitly set to zero by Fortran; other dE/dx values and errors are separate |
| `VdHit_signalToNoise`, `VdUnHit_signalToNoise` | No information in this MC; basic hit data are not thereby invalid |
| `Part_hacParticleCode`, `Part_vtxCode`, `Event_fillNumber` | No useful classification or metadata information in this sample |
| `GenPart_vertex.fCoordinates.fT`, `GenPart_tau` | No usable time or lifetime information |

## Truth links and legacy output

- 63.25% of Photon candidates lacked a direct Sim link. A missing link is not a missed photon.
- Distinguish conflicting Gen origins from missing origins. Several Sim descendants may merge into one reconstructed object, so asymmetric links alone are not a defect.
- No checked positive index was out of range. The `job_0` subset of 22,500 events had no array-length, Gen–Sim reciprocity or ancestry-cycle issue. This tests saved lineage, not detector efficiency.
- Zmumu is missing `job_16`, run −3016, event 1271. No duplicate events or empty Gen collections were found in that sample; the cause of the loss is unresolved.
- The separate legacy `t` TTree has `pid=-999` placeholders, artificial `ntpc=0/7`, incorrect `dpp` and negative `chi2ndf`. The paired datasets and current validation do not use it.

See [data conventions](../data.md) for indices and PDG versus DELPHI mass codes.
The source comparison used the
[native writer at `706d1dc`](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp);
the original executable's build hash was not verified.
