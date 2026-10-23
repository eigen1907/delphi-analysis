# Classification branch review

Audit of converter `47bc812` and the Jongwon dataset. **Pass means correct porting,
not validated physical quality.** Existing Florian files do not inherit later fixes.
Unverified means the external feature extraction code was unavailable.

The confirmed Raw MUID errors are separate from source-data outliers and the
`Part_charge` boundary case below. This review did not change the converter.

## Raw: `raw_sdst/Events`

| Branch | Status | Reason |
| --- | --- | --- |
| `TracRaw_theta` | Pass | Source angle (rad) and offset match. |
| `TracRaw_phi` | Pass | Source angle (rad) and offset match. |
| `TracRaw_invR` | Pass | Source curvature (1/cm) and offset match. |
| `TracRaw_charge` | Pass | Charge codes map to ±1; unknown maps to 0. |
| `Event_bFieldGevCm` | Pass | Fortran call and unit ratio match. |
| `TracRaw_impactRPhi` | Pass | Source impact parameter offset matches. |
| `TracRaw_impactZ` | Pass | Source impact parameter offset matches. |
| `TracRaw_trackLength` | Pass | Stores the absolute length as defined. |
| `TracRaw_chi2VD` | Pass | Source offset matches. |
| `TracRaw_ndfVD` | Pass | Source offset matches; out-of-range values become 0. |
| `TracRaw_chi2NoVD` | Pass | Copy is correct; negative/extreme source values remain unexplained. |
| `TracRaw_firstPointR` | Pass | First-point radius is calculated correctly. |
| `TracRaw_firstPointZ` | Pass | First-point Z is copied correctly. |
| `MtpcRaw_dEdx80Max` | Pass | MTPC source offset matches. |
| `MtpcRaw_dEdx80Sigma` | Pass | Stores the distribution width, not the measurement error. |
| `MtpcRaw_dEdx80Integrated` | Pass | MTPC source offset matches. |
| `MtpcRaw_dEdx65Max` | Pass | MTPC source offset matches. |
| `MtpcRaw_dEdx65Sigma` | Pass | Stores the distribution width, not the measurement error. |
| `MtpcRaw_zFitChi2` | Pass | Source offset matches. |
| `EmShower_energy` | Pass | EM shower energy is copied correctly. |
| `EmShower_position` | Pass | XYZ position (cm) is copied correctly. |
| `EmShower_nLayers` | Pass | HPC value is correct; FEMC has no layer information and uses −1. |
| `EmShower_detector` | Pass | HPC/FEMC detector code is correct. |
| `EmLayer_emShowerIdx` | Pass | Zero-based shower link is correct. |
| `HadShower_energy` | Pass | Hadronic shower energy is copied correctly. |
| `HadShower_direction` | Pass | Copy is correct; the value is an XYZ position (cm). |
| `HadShower_nHits` | Pass | Counts layer rows, not the total towers. |
| `HadHit_hadShowerIdx` | Pass | Zero-based shower link is correct. |
| `Stic_energyFromMain` | Pass | MAIN energy is copied correctly. |
| `Stic_directionFromMain` | Pass | Copy is correct; the value is momentum (GeV). |
| `Stic_numHitTowers` | Pass | Tower count decoding matches the source. |
| `VdAssocHit_tracRawIdx` | Pass | Raw track index link is correct. |
| `MuidRaw_tag` | Fail | Stores module number 24 instead of the muon tag. |
| `MuidRaw_looseChi2` | Fail | Stores the block bit pattern instead of χ². |
| `MuidRaw_hitPattern` | Fail | Stores a different pattern from before inefficiency correction. |
| `ElidRaw_tag` | Pass | Stored tag is correct; rows without an ELID module are omitted. |
| `ElidRaw_gammaConvTag` | Pass | Stored tag is correct; rows without an ELID module are omitted. |
| `HaidRaw_kaonDedx` | Pass | Decoded HAID offset and type match. |
| `HaidRaw_protonDedx` | Pass | Decoded HAID offset and type match. |
| `HaidRaw_kaonRich` | Pass | Decoded HAID offset and type match. |
| `HaidRaw_protonRich` | Pass | Decoded HAID offset and type match. |
| `HaidRaw_pionRich` | Pass | Decoded HAID offset and type match. |
| `HaidRaw_kaonCombined` | Pass | Decoded combined PID is copied correctly. |
| `HaidRaw_protonCombined` | Pass | Decoded combined PID is copied correctly. |
| `HaidRaw_richQuality` | Pass | Current int32 preserves the source status word. |
| `HaidRaw_thetaGas` | Pass | Corrected gas COMMON offset and length match. |
| `HaidRaw_thetaLiq` | Pass | Corrected liquid COMMON offset and length match. |
| `Vtx_position` | Pass | Stores XYZ from the original PV bank. |
| `Vtx_chi2` | Pass | Stores χ² from the original PV bank. |
| `Vtx_nOutgoing` | Pass | Stores the outgoing count from the original PV bank. |
| `nTracRaw` | Pass | Matches the stored Raw track count. |
| `nEmShower` | Pass | Matches the stored EM shower count. |
| `nHadShower` | Pass | Matches the stored hadronic shower count. |
| `nMtpcRaw` | Pass | Counts track rows with an MTPC module. |
| `nVdAssocHit` | Pass | Matches the stored associated VD hit count. |
| `nVtx` | Pass | Matches the stored Raw PV count. |
| `nStic` | Pass | Matches the stored STIC row count. |

## Native: `analysisLevel/Events`

The unknown charge code `999` was not observed in the checked sample.

| Branch | Status | Reason |
| --- | --- | --- |
| `Part_fourMomentum.fCoordinates.fX` | Pass | VECP Px is copied correctly. |
| `Part_fourMomentum.fCoordinates.fY` | Pass | VECP Py is copied correctly. |
| `Part_fourMomentum.fCoordinates.fZ` | Pass | VECP Pz is copied correctly. |
| `Part_fourMomentum.fCoordinates.fT` | Pass | Energy copy is correct; the 10.6 TeV outlier remains unexplained. |
| `Part_charge` | Fail (boundary) | ±1/0 work; int8 cannot preserve the unknown code 999. |
| `Part_hpcTotalShowerEnergy` | Pass | Associated HPC shower energies are summed correctly. |
| `Part_hpcNumLayers` | Pass | Stores the layer count of the strongest HPC shower. |
| `Part_hacTotalShowerEnergy` | Pass | Associated HAC shower energies are summed correctly. |
| `Part_hacNumTowers` | Pass | Matches the source Fortran tower accumulation. |
| `Part_sticShowerEnergy` | Pass | Corrected STIC COMMON length and offset match. |
| `Part_sticNumTowers` | Pass | Corrected STIC COMMON length and offset match. |
| `Trac_length` | Pass | Preserves length and the sign encoding PV-fit participation. |
| `Trac_impParRPhi` | Pass | Uses the correct scalar offset. |
| `Trac_impParZ` | Pass | Uses the correct scalar offset. |
| `Trac_chi2VD` | Pass | Decoded COMMON offset matches. |
| `Trac_ndfVD` | Pass | Decoded COMMON offset matches. |
| `Trac_chi2NoVD` | Pass | Copy is correct; negative/extreme source values remain unexplained. |
| `Trac_firstPointR` | Pass | First-point radius is copied correctly. |
| `Trac_firstPointZ` | Pass | First-point Z is copied correctly. |
| `Dedx_value` | Pass | Corrected dE/dx in MIP units is stored correctly. |
| `Muid_tag` | Pass | Stores the decoded muon tag for the PA ID. |
| `Muid_looseChi2` | Pass | Stores the correct decoded χ². |
| `Elid_tag` | Pass | Stores the correct decoded electron tag. |
| `Elid_gammaConversion` | Pass | Stores the correct decoded conversion tag. |
| `Haid_kaonDedx` | Pass | Decoded HAID offset and type match. |
| `Haid_protonDedx` | Pass | Decoded HAID offset and type match. |
| `Haid_kaonRich` | Pass | Decoded HAID offset and type match. |
| `Haid_protonRich` | Pass | Decoded HAID offset and type match. |
| `Haid_pionRich` | Pass | Decoded HAID offset and type match. |
| `Haid_richQuality` | Pass | Preserves the source quality word. |
| `Rich_theg` | Pass | Corrected gas COMMON offset and length match. |
| `Rich_thel` | Pass | Corrected liquid COMMON offset and length match. |
| `Haidc_pionTag` | Pass | KHAIDC offset and type match. |
| `Haidc_kaonTag` | Pass | KHAIDC offset and type match. |
| `Haidc_protonTag` | Pass | KHAIDC offset and type match. |
| `Haidc_heavyTag` | Pass | KHAIDC offset and type match. |
| `nPart` | Pass | Matches the native Part array length. |
| `nPhoton` | Pass | Counts neutral EM candidates; does not add PhotonConv. |
| `nV0` | Pass | Matches the decoded reconstructed V0 count. |
| `Event_chargedMult` | Pass | Original PA charged multiplicity. |
| `Event_neutralMult` | Pass | Original PA neutral multiplicity. |
| `Event_totalChargedEnergy` | Pass | Charged energy passing the Team4 selection. |
| `Event_totalEMEnergy` | Pass | Neutral EM energy. |
| `Event_totalHadronicEnergy` | Pass | Neutral hadronic energy. |
| `nVtx` | Pass | Final SKELANA vertex count. |
| `Vtx_chi2` | Pass | Final vertex χ²; its definition can differ from Raw PV. |
| `Vtx_nOutgoing` | Pass | Final vertex outgoing count; can differ from Raw. |
| `Event_hadronTagT4` | Pass | Source Team4 hadron tag conversion is correct. |

## External model features

Repeated branch names inherit the source-branch verdict. Selection, joins,
and calculations for `lead_*`, `sub_*`, and other derived features remain unverified.

| Branch | Status | Reason |
| --- | --- | --- |
| `nPart` | Pass | Source Part count is correct. |
| `n_charged` | Unverified | Charged selection and counting code unavailable. |
| `n_neutral` | Unverified | Neutral selection and counting code unavailable. |
| `nPhoton` | Pass | Source neutral EM candidate count is correct. |
| `nV0` | Pass | Source V0 count is correct. |
| `Event_chargedMult` | Pass | Source PA charged multiplicity is correct. |
| `Event_neutralMult` | Pass | Source PA neutral multiplicity is correct. |
| `Event_totalChargedEnergy` | Pass | Source selected charged energy is correct. |
| `Event_totalEMEnergy` | Pass | Source neutral EM energy is correct. |
| `Event_totalHadronicEnergy` | Pass | Source neutral hadronic energy is correct. |
| `Event_hadronTagT4` | Pass | Source Team4 tag is correct. |
| `n_extra` | Unverified | Definition and calculation of extra objects unavailable. |
| `photon_E_sum` | Unverified | Summed collection and selection are unverified. |
| `max_Muid_tag` | Unverified | Aggregation code unavailable; maximum of a bitmask needs interpretation. |
| `n_muid` | Unverified | Muon selection unavailable; array length is not the muon count. |
| `lead_E` | Unverified | Leading selection and energy handling are unverified. |
| `lead_charge` | Unverified | Leading join and unknown-charge handling are unverified. |
| `lead_hpcTotalShowerEnergy` | Unverified | Source branch passes; leading join is unverified. |
| `lead_hacTotalShowerEnergy` | Unverified | Source branch passes; leading join is unverified. |
| `lead_Dedx_value` | Unverified | Source branch passes; leading join is unverified. |
| `lead_Muid_tag` | Unverified | Native tag passes; leading join is unverified. |
| `lead_Elid_tag` | Unverified | Source tag passes; leading join is unverified. |
| `lead_Haidc_pionTag` | Unverified | Source tag passes; leading join is unverified. |
| `lead_Haidc_heavyTag` | Unverified | Source tag passes; leading join is unverified. |
| `lead_Haidc_kaonTag` | Unverified | Source tag passes; leading join is unverified. |
| `lead_Haid_kaonDedx` | Unverified | Source tag passes; leading join is unverified. |
| `lead_Haid_kaonCombined` | Unverified | Source branch exists and passes; leading join is unverified. |
| `lead_p` | Unverified | Leading selection and momentum calculation unavailable. |
| `lead_pt` | Unverified | Leading selection and transverse-momentum calculation unavailable. |
| `lead_theta` | Unverified | Leading selection and angle calculation unavailable. |
| `sub_E` | Unverified | Subleading selection and energy handling are unverified. |
| `sub_charge` | Unverified | Subleading join and unknown-charge handling are unverified. |
| `sub_hpcTotalShowerEnergy` | Unverified | Source branch passes; subleading join is unverified. |
| `sub_hacTotalShowerEnergy` | Unverified | Source branch passes; subleading join is unverified. |
| `sub_Dedx_value` | Unverified | Source branch passes; subleading join is unverified. |
| `sub_Muid_tag` | Unverified | Native tag passes; subleading join is unverified. |
| `sub_Elid_tag` | Unverified | Source tag passes; subleading join is unverified. |
| `sub_Haidc_pionTag` | Unverified | Source tag passes; subleading join is unverified. |
| `sub_Haidc_heavyTag` | Unverified | Source tag passes; subleading join is unverified. |
| `sub_Haidc_kaonTag` | Unverified | Source tag passes; subleading join is unverified. |
| `sub_Haid_kaonDedx` | Unverified | Source tag passes; subleading join is unverified. |
| `sub_Haid_kaonCombined` | Unverified | Source branch exists and passes; subleading join is unverified. |
| `sub_p` | Unverified | Subleading selection and momentum calculation unavailable. |
| `sub_pt` | Unverified | Subleading selection and transverse-momentum calculation unavailable. |
| `sub_theta` | Unverified | Subleading selection and angle calculation unavailable. |

## Evidence

Compared the C++ writer with `skelana.car` and `stdcdes.car` from the linked
DSTANA 161018 build, plus the SDST bank specification `shortdst.des`.

| Check | Jongwon sample |
| --- | --- |
| Part and event counts | 5,000 events per class, five classes |
| PID | 1,000 events per class, five classes |
| Calorimeters | All 60,000 ZKK training events |
| Tracks and vertices | 5,000 ZKK events |

Detailed source extracts, numerical results, and reproduction scripts are kept
locally in `.cache/porting-audit/` and are not tracked in Git.
See the [native](native.md) and [Raw](raw.md) branch references for the broader audit.
