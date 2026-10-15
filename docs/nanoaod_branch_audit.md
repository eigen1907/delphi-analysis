# Native branch 사전

대상: `20260828_florian` PHOTOS-FSR 원본 `nanoaod.root:Events`.
2026-10-13에 100파일·449,999이벤트를 검사했다. 샘플별 90,000개이며 Zmumu만 89,999개다.
346개 필드 중 수치 payload 310개를 검사한 원자료는 [CSV](nanoaod_numeric_scan.csv)에 보존한다.
NaN/Inf가 없어도 잘못된 word를 읽은 유한값은 남는다.

현재 dataset에서는 `analysisLevel/Events`로 읽는다.
아래 “수정”은 [converter 변경](gamma_rich_converter_fixes.md)을 뜻하며 **기존 Florian 파일에는 소급되지 않는다.**

## 주요 branch

| Branch / collection | Florian에서 확인한 내용 | 현재 사용·수정 상태 |
| --- | --- | --- |
| `GenPart_*` | 운동량, PDG, status, ancestry 확인. 시간/lifetime 없음 | 위치는 cm. Gen–Sim 연결 교정은 수정 문서 참고 |
| `Part_*`, `Trac_partIdx` | 저장 운동량과 Part 연결 사용 가능 | Native track 운동량은 `Trac → Part_fourMomentum`으로 읽음 |
| `Part_stic*`, `Photon_stic*` | COMMON stride 6/9 불일치; silicon은 잘못된 slot 6 | 수정 producer에서 stride 9, silicon slot 9로 교정 |
| `Photon_*` | Neutral EM 후보. STIC 에너지 오류가 후보 목록과 detector flag에도 영향 | 수정 후에도 Photon PID 확정 집합은 아님 |
| `PhotonConv_*` | 직접 Sim index 45,543개 모두 `-1`; 일부 daughter도 결측 | 유효 daughter/parent ancestry 사용. 수정 producer는 확인된 직접 Sim link만 채움 |
| `Rich_theg/sigg/nphg/nepg/flagg` | Gas COMMON 길이 오류. `Rich_sigg` nonzero 1,699,997개 모두 subnormal | 수정 producer는 decoded gas bank 사용 |
| `Rich_thel/sigl/nphl/nepl/flagl` | `nepl` nonzero 2,099,957개 모두 subnormal, `flagl` 타입 오류. 다른 liquid 값도 이전 이벤트 잔류 가능 | 수정 producer는 decoded liquid bank 사용 |
| `Haid_*`, `Haidc_*`, `Haide_*`, `Haidn_*`, `Haidr_*` | 별도 PID bank. 위 Rich 배열 오류를 모든 PID tag의 오류로 확대하지 않음 | Native `Haid_richQuality`는 원래 int32. `Haid_sign`은 charge가 아닌 PID method word |
| `Trac_perigee`, `Trac_weightMatrix` | 각각 `QTRAC(3..7)`, `(8..22)`를 읽음. 올바른 위치는 `(4..8)`, `(9..23)` | 이번 Gamma/RICH 수정 대상 아님. 이 필드로 운동량·오차를 복원하지 않음 |
| `Trac_impParRPhi`, `Trac_impParZ` | 올바른 slot 4/5 사용 | 위 perigee offset 오류와 구분 |
| `Trac_chi2NoVD`, `Trac_chi2VDHits` | 양수 ndf인데 음수 NoVD χ²인 경우 91,446개; VD-hit χ²에 큰 tail | Raw에도 나타남. 실패-fit sentinel이나 정상 χ²라고 일괄 해석하지 않음 |
| `SimVtx_position`, `SimVtx_masscode` | 좌표 대신 index 비트, masscode 대신 다른 slot 저장 | 수정 producer에서 `QVTX(6..8)`, `KVTX(5)`로 교정. 기존 integer ancestry는 별개 |
| `Event_totalChargedEnergy`, `Event_totalEMEnergy` | 각각 선택된 charged / neutral EM energy | 모든 Part 에너지 합과 같을 필요 없음 |
| `Elid`, `Muid`, `Btag`, `Uter`, `V0`, `Vtx` | 이 감사만으로 PID·fit 성능 전체를 검증한 것은 아님 | Branch별 정의와 유효성을 확인해 사용 |

## 빈 정보와 초기값

아래는 **Florian에서 관측한 값**이다. All-zero 자체를 producer 오류로 단정하지 않는다.

| Branch | 확인된 제한 |
| --- | --- |
| `BeamSpot_errorFlag`, `BeamSpot_x/y/z`, `BeamSpot_sigmaX/Y/Z` | 이벤트별 갱신 없이 초기값 0. Error flag 0도 유효한 beam spot을 보장하지 않음 |
| `V0Hyp_kind/probability/vertex/momentum`, `V0Hyp_impactXY/impactXYErr/impactZ/impactZErr/impact3D/impact3DErr` | `nV0Hyp=8*nV0`이지만 물리량 모두 0. SDST version 105는 old hypothesis를 채우지 않음 |
| `Trac_resRPhiFirstPoint`, `Trac_errorResRPhiFirstPoint`, `Trac_resZFirstPoint`, `Trac_errorResZFirstPoint` | 잔차 정보 없음 |
| `Dedx_width`, `Dedx_gapWires` | Fortran에서 0으로 설정. 다른 dE/dx 값·오차와 구분 |
| `VdHit_signalToNoise`, `VdUnHit_signalToNoise` | 이번 MC에서 정보 없음. 기본 hit 정보까지 무효라는 뜻은 아님 |
| `Part_hacParticleCode`, `Part_vtxCode`, `Event_fillNumber` | 이번 데이터에서 분류·metadata 정보 없음 |
| `GenPart_vertex.fCoordinates.fT`, `GenPart_tau` | 시간/lifetime 분석용 정보 없음 |

## Truth와 legacy tree

- 직접 Photon→Sim link가 없는 후보는 63.25%였다. 직접 link 부재를 photon 미검출로 해석하지 않는다.
- 양방향 Gen origin의 conflict/no-origin을 구분한다. Sim descendants가 한 Reco에 합쳐질 수 있어 비대칭 자체는 오류가 아니다.
- 검사한 index에 양수 범위 초과는 없었다. `job_0` 22,500이벤트에서는 배열 길이·Gen–Sim reciprocity·cycle 문제가 없었다. 이는 saved-lineage 검사이며 검출기 효율 검증이 아니다.
- Index 기준과 PDG/DELPHI code 구분은 [사용 규칙](branch_usage_summary.md)을 따른다.
- Zmumu의 누락은 `job_16`, run −3016, event 1271이다. 같은 샘플의 event 중복·빈 Gen은 없었고 누락 원인은 미확정이다.

원본의 legacy `t` TTree는 `Events`와 별개다. `pid=-999` placeholder,
`ntpc=0/7` 가짜 count, 잘못 계산한 `dpp`, 음수 `chi2ndf`가 확인됐다.
현재 paired dataset/검증 코드는 이 tree를 입력으로 쓰지 않는다.

근거: [수치 CSV](nanoaod_numeric_scan.csv),
[당시 native writer `706d1dc`](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp),
[현재 수정 목록](gamma_rich_converter_fixes.md). 원본 실행 바이너리의 build hash까지 확인한 감사는 아니다.
