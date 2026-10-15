# Gamma/RICH converter 수정

범위: `706d1dc..47bc812`의 `ff6051a`(Gamma), `47bc812`(RICH).
DELSIM/DELANA와 원본 이벤트는 변경하지 않고 ROOT converter를 수정했다.
수정은 새로 변환한 파일에 적용된다. 기존 `20260828_florian` ROOT에는 소급되지 않는다.
원본 tree `Events`는 dataset에서 아래 두 tree로 저장된다.

## Native: `analysisLevel/Events`

| ROOT branch | 수정 전 → 수정 후 |
| --- | --- |
| `Part_paIdx` | **신규.** Part/VECP row → 실제 0-based PA ID, 없으면 `-1` |
| `Rich_theg`, `Rich_sigg`, `Rich_nphg`, `Rich_nepg`, `Rich_flagg`, `Rich_thel`, `Rich_sigl`, `Rich_nphl`, `Rich_nepl`, `Rich_flagl` | 잘못된 COMMON/auxiliary 접근 → 현재 이벤트의 decoded gas/liquid bank. Expected photon count 소수 보존 |
| `Part_sticShowerEnergy`, `Part_sticShowerTheta`, `Part_sticShowerPhi`, `Part_sticNumTowers`, `Part_sticChargedTag`, `Part_sticSiliconVertexPos` | STIC stride `6 → 9`; silicon field `KSTIC(6) → KSTIC(9)` |
| `Photon_sticShowerEnergy`, `Photon_sticNumTowers`, `Photon_sticChargedTag` | 같은 STIC stride 교정 |
| `SimPart_genIdx`, `GenPart_simIdx` | SH 순서 `ISTSH/ISHST` → Gen/LU 순서 `ISTLU/ILUST` |
| `SimVtx_masscode`, `SimVtx_position` | `KVTX(4) → (5)`, `QVTX(1..3) → (6..8)` |
| `PhotonConv_simPhotonIdx` | 잘못된 `KPHOC(3)-1` → 저장 conversion vertex의 reciprocal parent Part와 직접 Sim gamma 연결을 확인. 없으면 `-1` |

STIC 교정은 후보/flag에도 간접 영향을 줄 수 있으나 해당 선택식 자체는 바꾸지 않았다.
Native에 새로 추가된 이름은 `Part_paIdx` 하나다. 기존 branch 이름 삭제·변경은 없다.

## Raw SDST: `raw_sdst/Events`

| ROOT branch | 수정 전 → 수정 후 |
| --- | --- |
| `HaidRaw_thetaGas`, `HaidRaw_sigmaGas`, `HaidRaw_nphGas`, `HaidRaw_nepGas`, `HaidRaw_flagGas`, `HaidRaw_thetaLiq`, `HaidRaw_sigmaLiq`, `HaidRaw_nphLiq`, `HaidRaw_nepLiq`, `HaidRaw_flagLiq` | Native와 같은 corrected decoded RICH bank. 잘못 읽은 값/잔류값 제거, expected count 소수 보존 |
| `HaidRaw_paIdx` | 1-based VECP row → 실제 0-based PA ID, 없으면 `-1` |
| `HaidRaw_richQuality` | `int8 → int32`, 전체 status word 보존 |
| `ElidRaw_tag`, `ElidRaw_gammaConvTag` | ELID header → PA identification word의 bits 7..9, 10..12 |
| `ElidRaw_refitMomentum` | 다른 ELID word → decoded `QELID(3..5)` |
| `EmShower_nLayers`, `EmShower_layerPattern` | FEMC record 밖 읽기 → 없는 FEMC 값은 `-1`; HPC는 유지 |

HaID의 기존 PID tag 7개(`kaonDedx`, `protonDedx`, `kaonRich`, `protonRich`, `pionRich`,
`kaonCombined`, `protonCombined`) 계산식은 유지했다. `HaidRaw_paIdx` 교정으로 올바른 track 연결이 가능해졌다.
`EmShower_nLayerEnergies`/`EmLayer`의 FEMC loop는 원래 skip되어 값이 바뀌지 않았다.

## 설명만 정정

| 소스 / branch | 실제 의미 |
| --- | --- |
| Native `Part_pdgId`, `SimPart_pdgId` | 표준 PDG가 아닌 DELPHI mass code: gamma=21, e−=2, e+=−2 |
| Native `Photon_emEnergy` | FEMC 에너지. HPC 에너지는 별도 저장 |
| Native `Photon_sticChargedTag` | 현재 STIC 0/1/2/3=no info/hard photon/soft photon/electron. Legacy 음수 코드도 병기 |
| Native `PhotonConv_pairType` | 현재 21–24 양수 코드, legacy는 음수 가능 |
| Native/raw `GenPart_vertex` | 이번 DST 위치는 cm. 시간/lifetime은 제공되지 않음 |
| Raw `Stic_directionFromMain` | 단위벡터가 아닌 MAIN momentum(GeV) |
| Raw `nHaidRaw` | Neutral placeholder를 포함한 VECP row 수 |

이 설명 정정 자체는 값을 바꾸지 않는다. Native `Haid_richQuality`는 기존 int32 그대로다.
ROOT 6.36 전후 RNTuple API 대응도 branch 의미·이름을 바꾸지 않는다.

## 검증과 남은 제한

기존 ZKK PHOTOS-FSR 1,000이벤트의 SDST/full FADANA를 재변환했다.
최종 코드로 각각 100이벤트를 다시 변환해 이전 수정본과 모든 컬럼의 값이 일치함을 확인했다.

| 검사 | 결과 |
| --- | --- |
| RICH native/raw 동일 PA 비교 | SDST 5,388 rows, full FADANA 5,090 rows의 11필드 모두 일치. Subnormal·중복 PA·track 연결 누락 0 |
| Gamma 표본 | Photon 1,280개, conversion 68개. 직접 Sim gamma 25개, unresolved 43개 |
| Sim vertex | 1,643개 좌표 finite, Gen–Sim 위치 5,282쌍 일치 |
| Raw ELID | 같은 PA의 native 1,908 rows와 tag/momentum 일치. 다른 raw 177컬럼 유지 |
| FEMC/HPC layer | FEMC 316개 metadata는 `-1`; HPC 2,469개와 layer 12,775개 값 유지 |

- PA 0은 유효하다. PA와 VECP/Part row를 같은 숫자로 연결하지 않는다.
- RICH는 gas COMMON 길이 `6 → 5`, NEP auxiliary `INTEGER`, flag auxiliary `REAL`로 맞추고 `PSBEG`가 초기화한 decoded bank를 읽는다. 추가 reset은 필요하지 않았다.
- RICH coverage는 **양수 normal angle이 저장된 track 비율**이다. Quality-bit cut/PID 성능 검증은 추가하지 않았다.
- Photon은 EM 후보다. STIC의 부족한 truth link나 conversion의 unresolved origin을 새로 만들어 채우지 않는다. 분석은 저장 parent/daughter ancestry를 사용한다.
- Full FADANA의 빈 raw EM/ELID는 detector 반응 부재를 뜻하지 않는다. Raw EM reader는 SDST 형식용이다.
- Electron refit momentum에 약 714 GeV tail이 있어 물리 사용 전 추가 검증이 필요하다.
- Native track offset, raw MUID, TPC packed 값 등 **이 두 커밋 밖의 문제**는 [native](nanoaod_branch_audit.md)/[raw](nanoaod_raw_sdst_branch_audit.md) 사전에 남겼다.

검증 자료는 `data/20261014_stic_validation/`에 보존했다.
현재 paired dataset으로 실행하는 명령은 [README](../README.md), 선택 정의는 [Photon fiducial](fiducial_regions.md)을 따른다.
