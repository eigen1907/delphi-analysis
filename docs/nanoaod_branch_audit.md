# `nanoaod.root` branch audit

확인일: 2026-10-13. 대상은 `20260828_florian`의 PHOTOS-FSR 5개 샘플이다.
분석 코드와 ROOT 파일은 수정하지 않았다.

간단한 선택 기준: [Collection 단위 사용 요약](branch_usage_summary.md).

## 1. 범위와 요약

검사 경로:

```text
data/20260828_florian/20260828_100kTest_{sample}_photosFSR/final_root/job_*/nanoaod.root
```

| Sample | Files | Events |
| --- | ---: | ---: |
| Zee | 20 | 90,000 |
| Zmumu | 20 | 89,999 |
| Ztautau | 20 | 90,000 |
| ZKK | 20 | 90,000 |
| Zpipi | 20 | 90,000 |
| **Total** | **100** | **449,999** |

주 대상은 `Events` RNTuple이며, 같은 파일의 legacy `t` TTree도 별도로 확인했다.
100개 파일의 스키마는 같았다. 컨테이너 필드를 포함한 이름은 346개이고,
실제 수치 payload 310개를 전체 이벤트에 대해 검사했다.
NaN/Inf는 없었지만, 정수의 비트를 float로 읽는 오류는 NaN 없이 작은 유한값으로 나타났다.

| 상태 | Branch / group | 판단 |
| --- | --- | --- |
| **저장 오류 확인** | `Part_stic*`, `Photon_stic*` | STIC COMMON 배열 길이와 일부 필드 위치 불일치 |
| **파생 선택 영향** | `Photon_isSTIC`, `Photon_isFEMC`, Photon 후보 목록 | 잘못 읽은 STIC 에너지로 후보와 flag 결정 |
| **저장 오류 확인** | gas `Rich_*`, `Rich_nepl`, `Rich_flagl` | COMMON 배열 길이 또는 타입 불일치 |
| **이전 이벤트 값 잔류 확인** | `Rich_*` convenience arrays | 현재 이벤트에서 채워지지 않은 row에 과거 값이 남음 |
| **저장 오류 확인** | `Trac_perigee`, `Trac_weightMatrix` | 시작 위치가 각각 한 word 앞섬 |
| **저장 오류 확인** | `SimVtx_position`, `SimVtx_masscode` | 좌표 대신 정수 index 비트, masscode 대신 다른 slot 저장 |
| **초기값 저장 확인** | `BeamSpot_*` | 이벤트마다 갱신하지 않고 전부 0 저장 |
| **미충전 확인** | `V0Hyp_*` 물리량 | 해당 SDST 버전에서는 hypothesis가 채워지지 않음 |
| **추가 조사 필요** | `Trac_chi2NoVD` | 양의 ndf에서도 음수 χ²; raw에도 같은 값 |
| **정보 부재 / 의미 주의** | `PhotonConv_simPhotonIdx`, truth links, 여러 all-zero 필드 | 결측과 비정상 저장을 구분해야 함 |

아래의 “확인”은 ROOT 수치와 producer 소스의 대응을 확인했다는 뜻이다.
각 sample 생성에 사용된 실행 바이너리의 build hash까지 추적한 것은 아니다.
전체 수치 결과는 [nanoaod_numeric_scan.csv](nanoaod_numeric_scan.csv)에 있다.
`values`는 펼친 숫자 개수이며, 행렬/벡터에서는 오브젝트 개수와 다르다.
`subnormal`은 0이 아닌 float의 절댓값이 해당 타입의 최소 normal 값보다 작은 경우다.

## 2. STIC: 기존 문제 재확인

관련 branch:

```text
Part_sticShowerEnergy, Part_sticShowerTheta, Part_sticShowerPhi
Part_sticNumTowers, Part_sticChargedTag, Part_sticSiliconVertexPos
Photon_sticShowerEnergy, Photon_sticNumTowers, Photon_sticChargedTag
Photon_isSTIC, Photon_isFEMC
```

Fortran의 `LENSTC=9`와 C++의 `LENSTC=6`이 다르다.
두 번째 Part부터 다른 word를 읽게 된다. 예를 들어 에너지 자리에 정수 tower/tag의
비트를 float로 읽으면 `1.4e-45` 수준의 값이 생긴다.
다른 Part의 정상적인 에너지를 읽는 경우도 있어 작은 에너지 컷으로 해결할 수 없다.
[C++ STIC 선언](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-analysis/include/skelana/pscemf_hpc_hac_stic.hpp),
[Fortran `PSCSTC` 정의](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car).

추가로 `Part_sticSiliconVertexPos`는 slot 6을 읽지만 silicon 정보는 slot 9다.
slot 6은 combined veto tag다. 배열 길이만 수정해도 이 필드는 남는 문제다.
`Photon_sticChargedTag`의 설명에 있는 음수 photon code도 현재 Fortran 정의와 맞지 않는다.
현재 slot 5는 0/1/2/3 = no information / hard gamma / soft gamma / electron이다.
[Writer의 STIC 저장 및 Photon 선택](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L411),
[Fortran STIC reader](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

전체 데이터에서:

- Reco Photon 573,796개 중 STIC 에너지 양수는 11,309개이며, 그중 3,659개가 subnormal이다.
- Photon STIC tower/tag가 1,000,000보다 큰 경우는 각각 1,539개다.
- Part STIC tower/tag/silicon 값이 1,000,000보다 큰 경우는 각각 6,633개다.

대표 예: `Zpipi/job_0`, ROOT entry 80, Photon row 0 → Part row 3.
Reco 방향은 θ=79.98°인데 저장된 STIC shower 방향은 θ=175.45°다.

`Photon_isSTIC`은 잘못 읽은 에너지의 `>0`으로 결정되고,
`Photon_isFEMC`는 STIC/HPC flag를 배제해 결정된다.
따라서 **STIC flag뿐 아니라 Photon 후보 목록과 FEMC 분류도 영향을 받는다**.
Photon four-momentum은 Part에서 복사하므로 모든 Photon 운동량이 같은 방식으로 깨졌다는 뜻은 아니다.
이 결과로 STIC 검출기 자체의 성능 저하를 결론낼 수 없다.

## 3. RICH: 타입 오류 외에 이전 이벤트 값도 남음

### 3.1 Gas의 다섯 필드

C++ gas bank 길이는 6, Fortran은 5다. bank 뒤의 convenience arrays가 한 배열씩 밀린다.

| 저장 branch | 실제로 읽는 Fortran 영역 |
| --- | --- |
| `Rich_theg` | `SIGG`: Cherenkov angle이 아니라 angular resolution |
| `Rich_sigg` | `NPHG`: 정수 photon count를 float 비트로 해석 |
| `Rich_nphg` | `NEPG`: observed count가 아니라 expected count |
| `Rich_nepg` | `FLAGG`: expected count가 아니라 real sector flag |
| `Rich_flagg` | `PSCGRC` COMMON 영역 밖 |

`Rich_sigg`의 0이 아닌 값 1,699,997개가 **모두 subnormal**이다.
다섯 필드 모두 현재 이름의 물리량으로 사용할 수 없다.
특히 현재 RICH 코드의 gas “saved-angle coverage”는 실제 angle 저장 여부를 검사하는 지표가 아니다.
[C++ gas 선언](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-analysis/include/skelana/pscgrc.hpp),
[Fortran `PSCGRC` 정의](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car),
[Writer의 RICH 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L1307).

### 3.2 Liquid의 타입

Liquid bank 길이는 맞지만 두 convenience array의 타입이 다르다.

- `Rich_nepl`: Fortran INTEGER를 C++ float로 읽음. 0이 아닌 2,099,957개가 모두 subnormal.
- `Rich_flagl`: Fortran REAL을 C++ integer로 읽음. 현재 데이터에서는 전부 0이라 잘못된 nonzero flag를 실증할 수는 없지만 타입 불일치는 확인됨.
- `Rich_thel`, `Rich_sigl`, `Rich_nphl`: 위치와 타입은 맞음. 다만 아래 잔류값 문제가 있음.

[C++ liquid 선언](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-analysis/include/skelana/psclrc.hpp),
[Fortran `PSCLRC` 정의](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car).

`Haidr_*`와 `Haid_*Rich` PID tag는 별도 bank에서 읽는다.
이 convenience array 오류 때문에 PID tag 자체도 잘못 저장됐다고 결론내리지는 않았다.

### 3.3 이전 이벤트의 RICH 값

Fortran 이벤트 초기화는 bank를 지우지만 위 convenience arrays를 초기화하지 않는다.
RICH reader는 charged Part만 처리하고, HAID/radiator 정보가 없으면 값 갱신 전에 반환한다.
Writer는 모든 Part row의 convenience arrays를 저장한다.
따라서 **현재 이벤트에 RICH 정보가 없는 charged Part에도 이전 이벤트 값이 남을 수 있다**.
[Fortran `PSBEG`, `PSHAD`, `PSHHAD`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

대표 예: `Zee/job_0`, ROOT entry 15, Part row 0.

- charged Part, `Haid_sign=Haid_richQuality=Haidr_selectionFlag=Dedx_value=0`.
- 운동량은 이전 이벤트와 다름.
- `Rich_thel=0.664330125`, `Rich_sigl=0.017584676`, `Rich_nphl=16`이 이전 entry 14의 동일 row와 정확히 같음.
- 현재 Track/RICH selection을 통과하므로 neutral Part 제외만으로 해결되지 않음.

각 sample의 `job_0` 총 22,500 이벤트에서, charged Part 중 quality=0이면서
양수 liquid angle과 angle/sigma/nph triplet이 이전 이벤트의 같은 row와 정확히 같은 경우는 2,413개였다.
전체에서는 neutral Part 695,391개에 양수 liquid angle이 있었다.
이 695,391개 전부가 잔류값이라고 개별 입증한 것은 아니다.

결론: **현재 liquid saved-angle coverage도 검증된 RICH coverage로 해석할 수 없다**.
raw SDST의 대응 HaidRaw 값도 같은 convenience API를 사용하므로 독립적인 정상 대조군이 아니다.
Producer의 타입/길이와 이벤트 초기화를 함께 수정하고 다시 생성하는 것이 우선이다.

## 4. Track: perigee / matrix 시작 위치 오류

| Branch | Writer에서 읽음 | Fortran에서 필요한 위치 |
| --- | --- | --- |
| `Trac_perigee` | `QTRAC(3..7)` | `QTRAC(4..8)` = d0, z0, θ, φ, 1/R |
| `Trac_weightMatrix` | `QTRAC(8..22)` | `QTRAC(9..23)` = packed weights 15개 |

Perigee의 첫 값은 integer Sim index의 비트를 읽고, curvature는 빠진다.
Matrix는 curvature로 시작하고 마지막 weight를 누락한다.
[Writer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L1040),
[Fortran `PSCTRA` 정의](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car).

대표 예: `ZKK/job_0`, ROOT entry 1, 첫 charged Part.

```text
native perigee: [0, -0.0941329, -0.389553, 1.434194, 1.157057]
raw perigee:    [-0.0941329, -0.389553, 1.434194, 1.157057, -8.02628e-5]
```

Native matrix 첫 값도 `-8.02628e-5`다. 원래 첫 weight는 `4495876.5`이며,
마지막 `351608832000`은 빠져 있다.

현재 native Track matching은 `Trac_partIdx → Part_fourMomentum`을 사용하므로 이 오류의 직접 대상은 아니다.
별도 `Trac_impParRPhi`/`Trac_impParZ`는 올바른 slot 4/5를 읽는다.
모든 `Trac_*`를 일괄적으로 잘못됐다고 분류하면 안 된다.

## 5. Sim vertex: 좌표가 정수 index의 비트

관련 branch는 `SimVtx_position.fCoordinates.fX/Y/Z`, `SimVtx_masscode`다.
Writer는 좌표를 `QVTX(1..3)`에서 읽지만, 실제 위치는 `QVTX(6..8)`이다.
Masscode도 `KVTX(4)`가 아니라 `KVTX(5)`다.
[Writer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L963),
[Fortran `PSCVTX` 정의](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car),
[Fortran `PSHSIM` 저장](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

전체 1,437,015 Sim vertices에서 float 좌표를 int32 비트로 읽으면 다음이 **모두 정확히 일치**했다.

```text
bits(X) = SimVtx_outgoingIdx + 1
bits(Y) = SimVtx_incomingIdx + 1
bits(Z) = SimVtx_nOut
```

`SimVtx_masscode`도 전부 0이다. Vertex 위치나 물질/변환 지점 분석에는 사용할 수 없다.
반면 해당 integer index들은 이 좌표 오류와 별개다.
현재 Gen–Sim–Reco ancestry traversal은 index를 사용하므로 좌표 오류 때문에 자동으로 무효가 되는 것은 아니다.

## 6. BeamSpot / V0Hyp: 초기값 또는 빈 payload

### BeamSpot

`BeamSpot_errorFlag`, `BeamSpot_x/y/z`, `BeamSpot_sigmaX/Y/Z`는
449,999 이벤트에서 모두 0이다.
`fillBeamSpot()`은 이벤트 처리에서 호출되지 않고 초기화 때 한 번 호출된다.
따라서 errorFlag=0도 유효한 이벤트 beam spot을 뜻하지 않는다.
Legacy `t.bx/by/ebx/eby`는 별도로 이벤트마다 저장하므로 같은 결론을 적용하지 않았다.
[Writer 초기화](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L123).

### V0Hyp

전체에서 `nV0Hyp=8*nV0`이며 1,461,968 hypothesis rows가 있지만,
다음 물리량은 모두 0이다.

```text
V0Hyp_kind, V0Hyp_probability, V0Hyp_vertex, V0Hyp_momentum
V0Hyp_impactXY, V0Hyp_impactXYErr, V0Hyp_impactZ, V0Hyp_impactZErr
V0Hyp_impact3D, V0Hyp_impact3DErr
```

현재 `Event_shortDstVersion=105`이고, Fortran은 해당 old hypothesis를 버전 ≤102에서만 채운다.
`v0Idx`/slot이 있어도 fit 정보가 있는 것은 아니다.
`kind != -1`만 요구하면 빈 kind=0 rows가 통과한다.
[Fortran `PSHRV0`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

## 7. 원인이 아직 확정되지 않은 Track χ²

Charged Part 1,417,964개 중 `Trac_chi2NoVD<0`은 652,654개다.
그중 **91,446개는 `Trac_ndfNoVD>0`인데도 χ²가 음수**다.
단순히 “fit 없음” sentinel만으로 설명되지 않는다.

예: `Zmumu/job_0`, ROOT entry 4에서 χ² `[-3.765205, -0.703876]`, ndf `[3, 1]`.
대응 raw SDST에도 같은 값이 있어 native writer만의 오류로 단정할 수 없다.
`Trac_chi2VD`에서는 charged Part의 음수 값이 없었다.
`Trac_chi2VDHits`에는 최대 `8.50e16`의 큰 tail도 있지만 원인은 확인하지 못했다.

NoVD fit 결과의 upstream 저장/정의를 확인하기 전에는 이 값으로 Track quality cut이나 ML feature를 만들지 않는 것이 좋다.

## 8. 정보가 없거나, 이름만 보고 해석하면 위험한 필드

### Photon conversion

- `PhotonConv_simPhotonIdx`: 45,543개 모두 -1. 이 직접 Sim link는 사용할 수 없다.
- 원본 `KPHOC(3)`는 Reco PA index를 찾는데 writer는 Sim index라고 설명한다. 현재 값이 모두 -1이므로 **positive index의 실제 오매핑은 이번 데이터에서 입증하지 못했다**.
- `PhotonConv_pairType`: 실제 +21/+23/+24. 최신 SDST의 양수 code와 맞으며, 음수 code라는 설명이 오래된 것이다.
- `PhotonConv_secondDaughterIdx=-1`: 9,619개. 한 TPC track conversion type에서는 예상 가능한 결측이다.
- `PhotonConv_firstDaughterIdx=-1`: 33개. 이 후보들은 daughter 기반 association이 제한되며 별도 확인이 필요하다.

[Conversion writer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp),
[Fortran `PSHPHC`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

### 전체가 0인 추가 필드

| Branch | 이번 데이터에서의 해석 |
| --- | --- |
| `Trac_resRPhiFirstPoint`, `Trac_errorResRPhiFirstPoint`, `Trac_resZFirstPoint`, `Trac_errorResZFirstPoint` | 잔차 정보 없음; 저장 원인은 추가 확인 필요 |
| `Dedx_width`, `Dedx_gapWires` | Fortran에서 0으로 설정; 측정값으로 사용 불가 |
| `VdHit_signalToNoise`, `VdUnHit_signalToNoise` | 이번 MC에서 정보 없음; 원인은 추가 확인 필요 |
| `Part_hacParticleCode`, `Part_vtxCode` | 이번 데이터에서 분류 정보 없음 |
| `GenPart_vertex.fCoordinates.fT`, `GenPart_tau` | 시간/lifetime 분석용 정보 없음 |
| `Event_fillNumber` | MC에서 fill metadata 없음 |

**All-zero라는 이유만으로 producer bug라고 판단하지 않았다.**
또한 `Photon_emEnergy`는 설명과 달리 HPC를 합친 EM 에너지가 아니라 FEMC/EMF의 `QEMF(8)`이다.
HPC 에너지는 `Photon_hpcTotalShowerEnergy`에 별도로 있다.
[Writer field 정의](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L453).

## 9. Truth link: 결측과 방향 불일치를 저장 오류와 구분

전체 Reco Photon의 직접 `Photon_partIdx → Part_simIdx` 결과:

| Direct Sim species | Count | 전체 Photon 대비 |
| --- | ---: | ---: |
| gamma (DELPHI code 21) | 183,169 | 31.92% |
| e− (2) | 3,274 | 0.57% |
| e+ (-2) | 3,361 | 0.59% |
| Other | 21,063 | 3.67% |
| No valid direct Sim link | 362,929 | 63.25% |
| **Total** | **573,796** | **100%** |

Ancestry resolver에서 forward/reverse Gen origin을 비교하면:

| Status | Count |
| --- | ---: |
| Agree | 198,319 |
| Forward only | 396 |
| Reverse only | 5,686 |
| Conflict | 12,153 |
| No origin | 357,242 |

Part→Sim과 Sim→Part가 상호 되돌아오지 않는 경우도 각각
60,644/1,563,174 valid links, 457,200/1,959,730 valid links였다.
여러 Sim descendants가 하나의 Reco에 연결될 수 있으므로 비대칭 자체를 오류로 볼 수 없다.
다만 conflict와 no-origin을 특정 Gen에 억지로 배정하면 안 된다.

전체 파일에서 검사한 `Part_simIdx`, `SimPart_partIdx`, `Photon_partIdx`에는 양수 범위 초과가 없었다.
각 `job_0` 총 22,500 이벤트의 정밀 검사에서는 검사 대상 배열 길이 불일치,
Gen↔Sim reciprocity 오류, ancestry cycle이 없었다.
이 subset의 stable Gen photons 89,763개 모두 GenPart_simIdx가 유효했다.
이는 **saved-lineage coverage**이고 검출기 transport/reconstruction 효율이 아니다.
Nested Gen anchor가 있으므로 기존 resolver처럼 **첫 유효 Gen anchor에서 멈춰야 한다**.

주의할 저장 규칙:

- `GenPart_pdgId`는 표준 PDG지만 `Part_pdgId`/`SimPart_pdgId`는 DELPHI mass code다. gamma=21, e±=±2.
- `Vtx_incomingIdx/outgoingIdx`는 1-based, 0 sentinel. 대응 `SimVtx_*`는 0-based, -1 sentinel.
- Native `Trac`/`Haid`/`Rich` rows는 neutral placeholder를 포함한 Part 순서다.
- `Haid_sign`은 charge가 아니라 PID method word다. `Haid_richQuality`의 큰 integer bit 값만으로 오류라고 볼 수 없다.
- `Event_totalChargedEnergy`는 선택된 charged energy, `Event_totalEMEnergy`는 neutral EM energy다. 모든 Part 에너지 합과 같아야 하는 것은 아니다.

## 10. 같은 파일의 legacy `t` TTree

아래는 `Events` branch와 별개로 확인한 문제다. 대상 particle rows는 1,797,052개다.

| `t` branch | 확인 결과 |
| --- | --- |
| `pid` | 전부 -999 placeholder; 실제 PID로 사용 불가 |
| `ntpc` | 실제 TPC hit 수가 아니라 0 또는 literal 7; 현재 charge를 대입하기 전 buffer로 결정 |
| `dpp` | curvature/weight slot을 잘못 읽어 계산; 정상 momentum resolution으로 사용 불가 |
| `chi2ndf` | -1 이외 음수 89,387개; NoVD χ² 문제를 이어받음 |

`ntpc`는 charged인데 0인 row가 175,371개, neutral인데 7인 row가 177,317개였다.
[Legacy `fillParticles()`](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-nanoaod/src/nanoaod_writer.cpp#L1523).

## 11. 현재 study에 미치는 영향 / 다음 확인

- **Gamma:** STIC/FEMC flag와 Photon 후보 선택은 producer 수정 후 재검사가 필요하다.
  Part/Gen 운동량과 integer ancestry를 쓰는 결과에는 같은 좌표/배열 오류가 직접 적용되지는 않지만,
  truth 결측과 conflict가 있어 truth-associated 효율을 검출기의 전체 효율로 단정할 수 없다.
- **Track:** Part 운동량 기반 matching과 별개로 perigee, weight matrix, NoVD χ²를 사용하는 분석은 재검토가 필요하다.
- **RICH:** 현재 gas/liquid 저장량과 coverage를 물리 성능으로 해석하지 말고, producer 수정 및 재생성이 필요하다.
- **Vertex / V0:** Sim 위치와 V0 hypothesis payload는 현재 파일로 물리 분석에 사용할 수 없다.

Zmumu의 1-event 부족은 `job_16`, run -3016의 event 1271이 없는 것으로 확인됐다.
같은 sample 안에 event ID 중복이나 Gen이 빈 이벤트는 없었다. 누락 원인은 branch 검사만으로 확정할 수 없다.

이번 감사는 native `nanoaod.root`가 우선이다. Raw SDST는 일부 원인 판별을 위한 대응값만 비교했다.
전체 raw schema/branch 건강 상태, NoVD χ²의 upstream 의미, 실제 생성 바이너리 provenance는 남은 확인 사항이다.
소스 기준은 delphi-nanoaod `706d1dc…`와 dstana `v270826`이다.

후속 검사: [nanoaod_raw_sdst.root branch audit](nanoaod_raw_sdst_branch_audit.md).
