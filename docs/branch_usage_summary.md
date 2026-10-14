# Branch 사용 요약 — collection 단위

기준: `20260828_florian` PHOTOS-FSR 5개 샘플, 2026-10-14까지의
[nanoaod 감사](nanoaod_branch_audit.md)와 [raw SDST 감사](nanoaod_raw_sdst_branch_audit.md).
`EmShower`는 `EmShower_*` 전체를 가리킨다.

| 판단 | 의미 |
| --- | --- |
| 사용 가능 | 아래에 적은 용도로 사용할 근거가 있음 |
| 부분 사용 | 사용할 수 있는 정보와 제외할 정보가 섞여 있음 |
| 사용 보류 | 현재 물리량/효율/PID 결과에 사용하기 어려움 |
| 사용 불가 | 해당 정보가 초기값 또는 빈 payload뿐임 |
| 추가 검증 | 특정 오류는 확인되지 않았지만 물리 의미/성능까지 검증하지 않음 |

“사용 가능”은 검출기 성능이나 PID 성능의 검증 완료를 뜻하지 않는다.

## `nanoaod.root`

| Collection | 판단 | 사용할 부분 / 제한 |
| --- | --- | --- |
| GenPart | 부분 사용 | 운동량, PDG, status, ancestry 사용. 시간/lifetime 제외 |
| Part | 부분 사용 | 운동량, charge, index 사용. STIC 정보 제외; truth 결측·충돌 구분 |
| SimPart | 사용 가능 | 저장된 운동량, species, ancestry 사용. Truth link의 결측·방향 불일치 구분 |
| SimVtx | 부분 사용 | Ancestry의 integer index만 사용. 좌표와 masscode 제외 |
| Trac | 부분 사용 | Part 연결과 직접 impact parameter 사용. Perigee, weight matrix, NoVD χ² 제외; VD-hit χ²도 추가 확인 |
| Photon | 부분 사용 | 저장 후보의 운동량·연결을 조사하는 용도. 후보 목록과 STIC/FEMC flag가 영향을 받아 최종 photon 효율 해석 보류 |
| PhotonConv | 부분 사용 | 운동량과 유효한 daughter 연결 사용. 직접 Sim link 없음; unresolved daughter는 별도 처리 |
| Rich | 사용 보류 | Gas/liquid 저장 오류와 이전 이벤트 값 잔류. Coverage와 PID response에 사용 보류 |
| BeamSpot | 사용 불가 | 이벤트별 갱신 없이 초기값만 저장 |
| V0Hyp | 사용 불가 | Hypothesis 물리량이 채워지지 않음 |
| Dedx | 부분 사용 | 값·오차 등은 결측 sentinel을 구분해 사용. Width/gap 정보 없음 |
| VdHit | 부분 사용 | 기본 hit 정보 탐색 가능. Signal-to-noise 없음 |
| VdUnHit | 부분 사용 | 기본 hit 정보 탐색 가능. Signal-to-noise 없음 |
| Event | 부분 사용 | Event ID와 정의를 확인한 요약량 사용. Fill 정보 없음; 에너지 요약은 전체 Part 합이 아님 |
| Haid, Haidc, Haide, Haidn, Haidr | 추가 검증 | 별도 PID bank. Rich 보조 배열 오류를 이 PID tags 전체의 오류로 확대하지 않음 |
| Elid, Muid, Btag, Uter, V0, Vtx | 추가 검증 | 이번 감사만으로 collection 전체의 PID/fit 성능을 검증하지 않음 |

Native Track 운동량은 `Trac → Part` 연결로 사용한다.
`SimVtx`의 좌표 오류는 integer index를 쓰는 ancestry traversal과 구분한다.
`Part`/`SimPart` species는 DELPHI mass code, `GenPart` species는 표준 PDG다.

## `nanoaod_raw_sdst.root`

| Collection | 판단 | 사용할 부분 / 제한 |
| --- | --- | --- |
| GenPart | 부분 사용 | 운동량, PDG, status, ancestry 사용. 이번 SDST vertex는 cm; 시간/lifetime 없음 |
| Stic | 사용 가능 | 에너지와 MAIN momentum 사용. Direction이라는 이름의 값은 unit vector가 아님 |
| EmShower | 부분 사용 | 에너지, 위치, detector ID 사용. FEMC layer metadata 제외; 전체가 photon PID 후보는 아님 |
| EmLayer | 부분 사용 | HPC layer/index 정보 사용. Layer 에너지 합을 보정된 shower 에너지와 동일시하지 않음 |
| TracRaw | 부분 사용 | Perigee와 weight matrix 사용. NoVD χ² 제외; VD-hit χ², 드문 invalid θ/zero curvature 구분 |
| MtpcRaw | 부분 사용 | dE/dx 저장값 탐색 가능. Packed counters는 현재 integer 해석으로 사용 불가; χ²도 추가 확인 |
| ElidRaw | 사용 보류 | Electron/conversion tag와 refit momentum이 잘못 저장됨. Bank 존재와 PA index는 별개 |
| MuidRaw | 사용 보류 | Muon tag, χ², hit pattern이 잘못 저장됨. Bank 존재와 PA index는 별개 |
| HaidRaw | 부분 사용 / RICH 보류 | RICH response·quality와 현재 Track join 사용 보류. 별도 dE/dx/PID tags까지 모두 오류라고 결론내리지는 않음 |
| TrackElement | 사용 불가 | 모든 record에서 빈 collection |
| VdAssocHit | 부분 사용 | 기본 hit/index 정보 탐색 가능. Signal-to-noise 없음 |
| VdUnassocHit | 부분 사용 | 기본 hit 정보 탐색 가능. Signal-to-noise 없음 |
| Event | 부분 사용 | Event ID, beam spot, simulated PV 사용. Setup records 구분; fill 정보 없음 |
| HadShower, HadHit, Vtx | 추가 검증 | 기본 저장 위치/일부 index는 확인. Collection 전체의 물리 의미·fit 성능은 추가 검증 |

Raw에는 native의 `Part`, `SimPart`, `SimVtx`, `Photon`, `PhotonConv`와
Gen–Sim truth-link 정보가 없다.
Raw calorimeter 후보와 Gen의 비교는 angular study로 할 수 있다.

**Raw PA index와 VECP/native Part index를 같은 숫자로 join하면 안 된다.**
현재 raw RICH Track association은 이 문제가 있어 물리 결과로 사용 보류다.

## 현재 분석에서 우선 사용할 조합

| Study | 우선 사용할 collection | 제한 |
| --- | --- | --- |
| Gen / ISR 정의 | 두 파일의 GenPart | PDG, status, ancestry 기반 |
| Native Track | Part + Trac 연결 | 잘못된 perigee/matrix와 의심되는 fit quality 제외 |
| Raw Track | TracRaw | 유효한 curvature/방향을 구분; NoVD χ² 제외 |
| Raw calorimeter Gamma | EmShower + Stic | Shower 후보의 angular study. Photon PID 효율과 구분 |
| Native conversion | PhotonConv + 유효한 daughter Part | 직접 Sim link 대신 daughter ancestry; 결측·모호성 구분 |
| Native truth ancestry | Part + SimPart + SimVtx index + GenPart | 첫 Gen anchor, 결측·충돌 처리 |
| RICH response / coverage | 사용 보류 | Producer 저장과 raw Track join 수정 후 재검증 |

Native 파일의 legacy `t` TTree는 별도 형식이다.
이 감사에서 확인된 PID placeholder, 가짜 TPC count, 잘못된 momentum error,
의심되는 χ²는 분석 입력으로 사용 보류다.
