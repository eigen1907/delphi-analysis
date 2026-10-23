# 분류 입력 branch 포팅 검토

현재 수정본(`47bc812`, Jongwon dataset) 기준이다.
**정상은 포팅이 맞다는 뜻**이며 원본 값의 물리적 품질까지 보장하지는 않는다.
계산 코드를 확인하지 못한 파생 feature는 `확인불가`로 표시했다.
Florian 기존 파일에는 현재 코드의 수정이 소급되지 않는다.

## Raw — `raw_sdst/Events`

| Branch 이름 | 정상/비정상 | 이유 |
| --- | --- | --- |
| `TracRaw_theta` | 정상 | 원본 각도(rad) 저장 위치 일치 |
| `TracRaw_phi` | 정상 | 원본 각도(rad) 저장 위치 일치 |
| `TracRaw_invR` | 정상 | 원본 곡률(1/cm) 저장 위치 일치 |
| `TracRaw_charge` | 정상 | ±1 변환 정상. Unknown은 0으로 저장 |
| `Event_bFieldGevCm` | 정상 | Fortran 호출·단위 비율 일치 |
| `TracRaw_impactRPhi` | 정상 | 원본 impact parameter 위치 일치 |
| `TracRaw_impactZ` | 정상 | 원본 impact parameter 위치 일치 |
| `TracRaw_trackLength` | 정상 | 정의대로 길이의 절댓값 저장 |
| `TracRaw_chi2VD` | 정상 | 원본 저장 위치 일치 |
| `TracRaw_ndfVD` | 정상 | 저장 위치 일치. 범위 밖 값은 0 처리 |
| `TracRaw_chi2NoVD` | 정상 | 복사는 정상이나 음수·극단값 원인 미확정 |
| `TracRaw_firstPointR` | 정상 | 첫 측정점의 반경 계산 정상 |
| `TracRaw_firstPointZ` | 정상 | 첫 측정점 Z 저장 정상 |
| `MtpcRaw_dEdx80Max` | 정상 | MTPC 원본 위치 일치 |
| `MtpcRaw_dEdx80Sigma` | 정상 | 분포 폭 저장 정상. 측정 오차와는 다름 |
| `MtpcRaw_dEdx80Integrated` | 정상 | MTPC 원본 위치 일치 |
| `MtpcRaw_dEdx65Max` | 정상 | MTPC 원본 위치 일치 |
| `MtpcRaw_dEdx65Sigma` | 정상 | 분포 폭 저장 정상. 측정 오차와는 다름 |
| `MtpcRaw_zFitChi2` | 정상 | 원본 저장 위치 일치 |
| `EmShower_energy` | 정상 | EM shower 에너지 복사 정상 |
| `EmShower_position` | 정상 | XYZ 위치(cm) 복사 정상 |
| `EmShower_nLayers` | 정상 | HPC 값 정상. FEMC는 정보가 없어 −1 |
| `EmShower_detector` | 정상 | HPC/FEMC detector code 정상 |
| `EmLayer_emShowerIdx` | 정상 | 0-based shower 연결 정상 |
| `HadShower_energy` | 정상 | Hadronic shower 에너지 복사 정상 |
| `HadShower_direction` | 정상 | 복사 정상. 실제 의미는 XYZ 위치(cm) |
| `HadShower_nHits` | 정상 | Layer row 수 저장 정상. Tower 총수와 다름 |
| `HadHit_hadShowerIdx` | 정상 | 0-based shower 연결 정상 |
| `Stic_energyFromMain` | 정상 | MAIN 에너지 복사 정상 |
| `Stic_directionFromMain` | 정상 | 복사 정상. 실제 의미는 운동량(GeV) |
| `Stic_numHitTowers` | 정상 | 원본 tower 수 decoding과 일치 |
| `VdAssocHit_tracRawIdx` | 정상 | Raw track index 연결 정상 |
| `MuidRaw_tag` | 비정상 | 실제 tag 대신 모듈 번호 24 저장 |
| `MuidRaw_looseChi2` | 비정상 | χ² 대신 block bit pattern 저장 |
| `MuidRaw_hitPattern` | 비정상 | 비효율 보정 전의 다른 pattern 저장 |
| `ElidRaw_tag` | 정상 | 저장된 tag 정상. ELID 모듈 없는 row는 누락 |
| `ElidRaw_gammaConvTag` | 정상 | 저장된 tag 정상. ELID 모듈 없는 row는 누락 |
| `HaidRaw_kaonDedx` | 정상 | Decoded HAID 위치·타입 일치 |
| `HaidRaw_protonDedx` | 정상 | Decoded HAID 위치·타입 일치 |
| `HaidRaw_kaonRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `HaidRaw_protonRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `HaidRaw_pionRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `HaidRaw_kaonCombined` | 정상 | Decoded combined PID 복사 정상 |
| `HaidRaw_protonCombined` | 정상 | Decoded combined PID 복사 정상 |
| `HaidRaw_richQuality` | 정상 | 현재 int32로 원본 status word 보존 |
| `HaidRaw_thetaGas` | 정상 | 수정된 gas COMMON 위치·길이 일치 |
| `HaidRaw_thetaLiq` | 정상 | 수정된 liquid COMMON 위치·길이 일치 |
| `Vtx_position` | 정상 | 원래 PV bank의 XYZ 저장 |
| `Vtx_chi2` | 정상 | 원래 PV bank의 χ² 저장 |
| `Vtx_nOutgoing` | 정상 | 원래 PV bank의 outgoing 수 저장 |
| `nTracRaw` | 정상 | 저장한 Raw track 수와 일치 |
| `nEmShower` | 정상 | 저장한 EM shower 수와 일치 |
| `nHadShower` | 정상 | 저장한 hadronic shower 수와 일치 |
| `nMtpcRaw` | 정상 | MTPC 모듈이 있는 track row 수와 일치 |
| `nVdAssocHit` | 정상 | 저장한 associated VD hit 수와 일치 |
| `nVtx` | 정상 | 저장한 Raw PV 수와 일치 |
| `nStic` | 정상 | 저장한 STIC row 수와 일치 |

## Native — `analysisLevel/Events`

| Branch 이름 | 정상/비정상 | 이유 |
| --- | --- | --- |
| `Part_fourMomentum.fCoordinates.fX` | 정상 | VECP의 Px 복사 정상 |
| `Part_fourMomentum.fCoordinates.fY` | 정상 | VECP의 Py 복사 정상 |
| `Part_fourMomentum.fCoordinates.fZ` | 정상 | VECP의 Pz 복사 정상 |
| `Part_fourMomentum.fCoordinates.fT` | 정상 | 에너지 복사 정상. 최대 10.6 TeV outlier 원인 미확정 |
| `Part_charge` | 비정상(경계값) | ±1/0은 정상. Unknown 999는 int8로 보존 불가 |
| `Part_hpcTotalShowerEnergy` | 정상 | Associated HPC shower 에너지 합 정상 |
| `Part_hpcNumLayers` | 정상 | 가장 강한 HPC shower의 layer 수 저장 |
| `Part_hacTotalShowerEnergy` | 정상 | Associated HAC shower 에너지 합 정상 |
| `Part_hacNumTowers` | 정상 | 원본 Fortran의 tower 누적 방식과 일치 |
| `Part_sticShowerEnergy` | 정상 | 수정된 STIC COMMON 길이·위치 일치 |
| `Part_sticNumTowers` | 정상 | 수정된 STIC COMMON 길이·위치 일치 |
| `Trac_length` | 정상 | 길이와 PV fit 참여 여부를 나타내는 부호 보존 |
| `Trac_impParRPhi` | 정상 | 올바른 scalar 위치 사용 |
| `Trac_impParZ` | 정상 | 올바른 scalar 위치 사용 |
| `Trac_chi2VD` | 정상 | Decoded COMMON 위치 일치 |
| `Trac_ndfVD` | 정상 | Decoded COMMON 위치 일치 |
| `Trac_chi2NoVD` | 정상 | 복사는 정상이나 음수·극단값 원인 미확정 |
| `Trac_firstPointR` | 정상 | 첫 측정점 반경 복사 정상 |
| `Trac_firstPointZ` | 정상 | 첫 측정점 Z 복사 정상 |
| `Dedx_value` | 정상 | 보정된 MIP 기준 dE/dx 저장 정상 |
| `Muid_tag` | 정상 | PA ID에서 decoded muon tag 저장 |
| `Muid_looseChi2` | 정상 | 올바른 decoded χ² 저장 |
| `Elid_tag` | 정상 | 올바른 decoded electron tag 저장 |
| `Elid_gammaConversion` | 정상 | 올바른 decoded conversion tag 저장 |
| `Haid_kaonDedx` | 정상 | Decoded HAID 위치·타입 일치 |
| `Haid_protonDedx` | 정상 | Decoded HAID 위치·타입 일치 |
| `Haid_kaonRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `Haid_protonRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `Haid_pionRich` | 정상 | Decoded HAID 위치·타입 일치 |
| `Haid_richQuality` | 정상 | 원본 quality word 보존 |
| `Rich_theg` | 정상 | 수정된 gas COMMON 위치·길이 일치 |
| `Rich_thel` | 정상 | 수정된 liquid COMMON 위치·길이 일치 |
| `Haidc_pionTag` | 정상 | KHAIDC 위치·타입 일치 |
| `Haidc_kaonTag` | 정상 | KHAIDC 위치·타입 일치 |
| `Haidc_protonTag` | 정상 | KHAIDC 위치·타입 일치 |
| `Haidc_heavyTag` | 정상 | KHAIDC 위치·타입 일치 |
| `nPart` | 정상 | Native Part 배열 길이와 일치 |
| `nPhoton` | 정상 | Neutral EM 후보 수. PhotonConv 합산 수는 아님 |
| `nV0` | 정상 | Decoded reconstructed V0 수와 일치 |
| `Event_chargedMult` | 정상 | 원래 PA 기준 charged multiplicity |
| `Event_neutralMult` | 정상 | 원래 PA 기준 neutral multiplicity |
| `Event_totalChargedEnergy` | 정상 | Team4 선택을 통과한 charged 에너지 |
| `Event_totalEMEnergy` | 정상 | Neutral EM 에너지 |
| `Event_totalHadronicEnergy` | 정상 | Neutral hadronic 에너지 |
| `nVtx` | 정상 | 최종 SKELANA vertex 수 |
| `Vtx_chi2` | 정상 | 최종 vertex χ². Raw PV와 정의가 다를 수 있음 |
| `Vtx_nOutgoing` | 정상 | 최종 vertex outgoing 수. Raw와 다를 수 있음 |
| `Event_hadronTagT4` | 정상 | 원본 Team4 hadron tag 변환 정상 |

`Part_charge`의 unknown 999 경계값은 검사한 표본에서 관측되지 않았다.

## 세 번째 목록 — 모델 입력 feature

동일한 branch 이름은 원본 branch의 판정이다.
`lead_*`, `sub_*` 등의 선택·연결·계산 과정은 추출 코드가 없어 확인하지 못했다.

| Branch 이름 | 정상/비정상 | 이유 |
| --- | --- | --- |
| `nPart` | 정상 | 원본 Part 수 정상 |
| `n_charged` | 확인불가 | Charged 선택·집계 코드 없음 |
| `n_neutral` | 확인불가 | Neutral 선택·집계 코드 없음 |
| `nPhoton` | 정상 | 원본 neutral EM 후보 수 정상 |
| `nV0` | 정상 | 원본 V0 수 정상 |
| `Event_chargedMult` | 정상 | 원본 PA charged multiplicity 정상 |
| `Event_neutralMult` | 정상 | 원본 PA neutral multiplicity 정상 |
| `Event_totalChargedEnergy` | 정상 | 원본 선택된 charged 에너지 정상 |
| `Event_totalEMEnergy` | 정상 | 원본 neutral EM 에너지 정상 |
| `Event_totalHadronicEnergy` | 정상 | 원본 neutral hadronic 에너지 정상 |
| `Event_hadronTagT4` | 정상 | 원본 Team4 tag 정상 |
| `n_extra` | 확인불가 | Extra의 정의·계산 코드 없음 |
| `photon_E_sum` | 확인불가 | 합산 collection·선택 조건 미확인 |
| `max_Muid_tag` | 확인불가 | 집계 코드 없음. Bitmask의 최댓값 해석 필요 |
| `n_muid` | 확인불가 | Muon 판정 조건 없음. 배열 길이와 muon 수는 다름 |
| `lead_E` | 확인불가 | Leading 선택·에너지 처리 미확인 |
| `lead_charge` | 확인불가 | Leading 연결·unknown charge 처리 미확인 |
| `lead_hpcTotalShowerEnergy` | 확인불가 | 원본 branch 정상. Leading 연결 미확인 |
| `lead_hacTotalShowerEnergy` | 확인불가 | 원본 branch 정상. Leading 연결 미확인 |
| `lead_Dedx_value` | 확인불가 | 원본 branch 정상. Leading 연결 미확인 |
| `lead_Muid_tag` | 확인불가 | Native tag 정상. Leading 연결 미확인 |
| `lead_Elid_tag` | 확인불가 | 원본 tag 정상. Leading 연결 미확인 |
| `lead_Haidc_pionTag` | 확인불가 | 원본 tag 정상. Leading 연결 미확인 |
| `lead_Haidc_heavyTag` | 확인불가 | 원본 tag 정상. Leading 연결 미확인 |
| `lead_Haidc_kaonTag` | 확인불가 | 원본 tag 정상. Leading 연결 미확인 |
| `lead_Haid_kaonDedx` | 확인불가 | 원본 tag 정상. Leading 연결 미확인 |
| `lead_Haid_kaonCombined` | 확인불가 | 원본 branch 존재·복사 정상. Leading 연결 미확인 |
| `lead_p` | 확인불가 | Leading 선택·운동량 계산 코드 없음 |
| `lead_pt` | 확인불가 | Leading 선택·횡운동량 계산 코드 없음 |
| `lead_theta` | 확인불가 | Leading 선택·각도 계산 코드 없음 |
| `sub_E` | 확인불가 | Subleading 선택·에너지 처리 미확인 |
| `sub_charge` | 확인불가 | Subleading 연결·unknown charge 처리 미확인 |
| `sub_hpcTotalShowerEnergy` | 확인불가 | 원본 branch 정상. Subleading 연결 미확인 |
| `sub_hacTotalShowerEnergy` | 확인불가 | 원본 branch 정상. Subleading 연결 미확인 |
| `sub_Dedx_value` | 확인불가 | 원본 branch 정상. Subleading 연결 미확인 |
| `sub_Muid_tag` | 확인불가 | Native tag 정상. Subleading 연결 미확인 |
| `sub_Elid_tag` | 확인불가 | 원본 tag 정상. Subleading 연결 미확인 |
| `sub_Haidc_pionTag` | 확인불가 | 원본 tag 정상. Subleading 연결 미확인 |
| `sub_Haidc_heavyTag` | 확인불가 | 원본 tag 정상. Subleading 연결 미확인 |
| `sub_Haidc_kaonTag` | 확인불가 | 원본 tag 정상. Subleading 연결 미확인 |
| `sub_Haid_kaonDedx` | 확인불가 | 원본 tag 정상. Subleading 연결 미확인 |
| `sub_Haid_kaonCombined` | 확인불가 | 원본 branch 존재·복사 정상. Subleading 연결 미확인 |
| `sub_p` | 확인불가 | Subleading 선택·운동량 계산 코드 없음 |
| `sub_pt` | 확인불가 | Subleading 선택·횡운동량 계산 코드 없음 |
| `sub_theta` | 확인불가 | Subleading 선택·각도 계산 코드 없음 |

## 검토 근거

실제 빌드에 연결된 DSTANA 161018의 `skelana.car`, `stdcdes.car`와 SDST bank 명세 `shortdst.des`를 대조했다.
Jongwon 5개 class 각 5,000이벤트에서 Part·event count를, 각 1,000이벤트에서 PID를 확인했다.
Calorimeter는 ZKK train 전체 60,000이벤트, track·vertex는 ZKK 5,000이벤트를 검사했다.
상세 근거·수치 JSON·재현 스크립트는 [`.cache/porting-audit/`](../.cache/porting-audit/)에 보존했다.
