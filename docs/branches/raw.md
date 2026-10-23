# Raw SDST branches

This audit covers the original `20260828_florian` PHOTOS-FSR files,
`nanoaod_raw_sdst.root:Events`. The scan on 2026-10-14 checked 100 files and
450,099 records: 449,999 physical events and one setup record per file. It
examined 171 numeric payloads among 185 fields and 135 collection-leaf lengths;
the [scan CSV](../scans/raw.csv) preserves the results.

Paired datasets store this tree as `raw_sdst/Events` and omit setup records.
The [converter fixes](../converter.md) apply to newly converted files, **not to
the existing Florian files**. For model inputs, see the
[classification review](classification.md).

## Branch status

| Branch / collection | Florian finding | Current status / use |
| --- | --- | --- |
| `GenPart_*` | PDG, status and parent match native in every physical event; ancestry indices are in range | Positions are in cm; no time or lifetime information |
| `ElidRaw_tag`, `ElidRaw_gammaConvTag`, `ElidRaw_refitMomentum` | Read ELID headers or unrelated values, including tag 25 and conversion tag 9040 | Fixed to PA identification bits and decoded electron momentum |
| `MuidRaw_tag`, `MuidRaw_looseChi2`, `MuidRaw_hitPattern` | Tag is module ID 24; chi-squared is a block pattern; hit pattern reads the wrong field | Still outside the Gamma/RICH fixes. Bank presence and PA links do not validate the PID payload |
| `MtpcRaw_packedPadsSectors`, `MtpcRaw_packedWiresHits` | Store integer bit patterns of packed REAL values; e.g. integer 1,195,722,496 represents float 50,503 | Still outside those fixes. Do not decode the stored integer directly with `%100`; other dE/dx measurements are separate |
| `EmShower_nLayers`, `EmShower_layerPattern` | Read beyond the four-word FEMC record | Fixed: unavailable FEMC values are `-1`; HPC values unchanged |
| `EmLayer_*` | HPC layer counts, indices and bit patterns agree | A difference between summed layer energy and corrected shower energy is not by itself a storage error |
| `HaidRaw_thetaGas/sigmaGas/nphGas/nepGas/flagGas` | Wrong gas length; all 1,700,905 nonzero `sigmaGas` values were subnormal | Fixed converter reads the decoded gas bank |
| `HaidRaw_thetaLiq/sigmaLiq/nphLiq/nepLiq/flagLiq` | All 2,103,405 nonzero `nepLiq` values were subnormal; wrong flag type and stale values | Fixed converter reads the decoded liquid bank |
| `HaidRaw_richQuality` | Narrowed to int8; e.g. native 256 became raw 0 | Fixed to int32. Other writer differences mean truncation alone does not explain every historical mismatch |
| `HaidRaw_paIdx`, `nHaidRaw` | Used 1-based VECP rows, including neutral placeholders | Fixed index is the actual 0-based PA ID or `-1`. Equal row counts do not make PA and Part indices interchangeable |
| `Stic_*` | MAIN energy/momentum offsets are correct; all 25,864 rows have positive finite energy matching momentum magnitude | `Stic_directionFromMain` is momentum in GeV. This is not a detector-efficiency test |
| `TracRaw_*` | Correct positions for five perigee and 15 weight values | Does not share the native composite-offset bug; rare invalid angles/curvatures and fit-quality meaning need separate checks |
| `TracRaw_chi2NoVD`, `TracRaw_chi2VDHits`, `MtpcRaw_zFitChi2` | 100,041 negative NoVD chi-squared values with positive ndf; large VD-hit tails; three negative z-fit values | Do not treat all as failed-fit sentinels. Packing and upstream definitions need further study |
| `Event_beamSpot*`, `Event_simPV*` | Updated per event; errorFlag is 0 for physical events and −1 for setup | Distinct from the unfilled native `BeamSpot_*` fields |
| `HadShower`, `HadHit`, `Vtx` | Some storage positions and indices checked | This scan alone does not validate every quantity or fit |

## Missing information and event selection

| Branch | Florian limitation |
| --- | --- |
| `TrackElement_*`, `nTrackElement` | Empty in every record; missing-bank or lookup cause unresolved |
| `VdAssocHit_signalToNoise`, `VdUnassocHit_signalToNoise` | Always zero; basic hit/index data are separate |
| `GenPart_productionTime`, `GenPart_properLifetime` | Filled with zero by the SDST unpacker |
| `Event_fillNumber` | No MC fill metadata |

Setup records have event number 0 and `nGenPart=0`. Removing them is different
from dropping physical events with no reconstructed photons or tracks. The scan
also found two tracks with theta greater than pi and 30 with zero inverse radius.

## Associations and photon studies

- `EmLayer_emShowerIdx`, `HadHit_hadShowerIdx`, `VdAssocHit_tracRawIdx` and `MtpcRaw_tracRawIdx` were valid 0-based indices in all scanned files.
- One PA may have several EM showers. Do not deduplicate clusters solely because they share a PA.
- Of 1,273,959 EM showers, 722,179 were associated with a charged raw-track PA. The full EM collection is not a photon-PID-selected collection.
- Raw has no native `Part`, `SimPart`, `SimVtx`, `Photon`, `PhotonConv` or reconstructed-to-Gen truth links. Gen ancestry and shower angular coverage can still be studied.
- Before the fix, joining PA and Part rows by equal numbers produced charge mismatches. Neither that join nor angular matching is an internal truth link.
- RICH quality 0 does not establish missing information. Do not interpret Florian RICH response or coverage until it is reconverted.

See [data conventions](../data.md) for index rules. The source comparison used
the [raw writer at `706d1dc`](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp);
the original executable's build hash was not verified.
