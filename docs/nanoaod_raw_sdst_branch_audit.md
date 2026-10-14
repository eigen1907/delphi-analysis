# `nanoaod_raw_sdst.root` branch audit

확인일: 2026-10-14. 대상은 `20260828_florian`의 PHOTOS-FSR 5개 샘플이다.
[Native nanoaod 감사](nanoaod_branch_audit.md)에 이어 raw 파일의 수치, 배열 길이, index와 producer 저장 위치를 확인했다.
분석 코드와 ROOT 파일은 수정하지 않았다.

간단한 선택 기준: [Collection 단위 사용 요약](branch_usage_summary.md).

## 1. 범위와 요약

```text
data/20260828_florian/20260828_100kTest_{sample}_photosFSR/final_root/job_*/nanoaod_raw_sdst.root
```

| Sample | Files | ROOT records | Physics events | Setup records |
| --- | ---: | ---: | ---: | ---: |
| Zee | 20 | 90,020 | 90,000 | 20 |
| Zmumu | 20 | 90,019 | 89,999 | 20 |
| Ztautau | 20 | 90,020 | 90,000 | 20 |
| ZKK | 20 | 90,020 | 90,000 | 20 |
| Zpipi | 20 | 90,020 | 90,000 | 20 |
| **Total** | **100** | **450,099** | **449,999** | **100** |

각 파일에 `nGenPart=0`, event number=0인 setup record가 하나 있다.
이 100개는 물리 이벤트 수에 포함하지 않았다. 선택한 Gen/Reco가 0개인 물리 이벤트를 버린 것은 아니다.
물리 이벤트의 run/event key는 native 파일과 모두 대응하며, Zmumu의 1-event 부족도 동일하다.

100개 파일의 스키마는 같았다. 컨테이너를 포함한 이름 185개, 수치 payload 171개를 전체 records에 대해 검사했다.
NaN/Inf는 없었다. 135개 collection payload leaf의 길이도 매 파일에서 해당 `nCollection`과 일치했다.
전체 수치 결과는 [nanoaod_raw_sdst_numeric_scan.csv](nanoaod_raw_sdst_numeric_scan.csv)에 있다.
`values`는 펼친 숫자 개수이므로 행렬에서는 오브젝트 개수와 다르다.
`subnormal`은 0이 아닌 float의 절댓값이 해당 타입의 최소 normal 값보다 작은 경우다.

| 상태 | Branch / group | 판단 |
| --- | --- | --- |
| **저장 위치 오류 확인** | `ElidRaw_tag`, `ElidRaw_gammaConvTag`, `ElidRaw_refitMomentum` | ELID header/다른 물리량을 electron ID와 momentum으로 저장 |
| **저장 위치 오류 확인** | `MuidRaw_tag`, `MuidRaw_looseChi2`, `MuidRaw_hitPattern` | MUID header/다른 tag를 잘못 저장 |
| **타입 해석 오류 확인** | `MtpcRaw_packedPadsSectors`, `MtpcRaw_packedWiresHits` | real로 저장된 packed 숫자를 integer 비트로 읽음 |
| **저장 위치 오류 확인** | FEMC의 `EmShower_nLayers`, `EmShower_layerPattern` | 4-word shower record 바깥에서 읽음 |
| **Native와 공통 오류** | gas `HaidRaw_*`, liquid `nepLiq/flagLiq` | 같은 잘못된 RICH COMMON/API 사용 |
| **이전 이벤트 값 잔류** | RICH convenience arrays | native처럼 이벤트별 초기화가 불완전함 |
| **Raw에서 추가 정보 손실** | `HaidRaw_richQuality` | status word를 int8로 축소 |
| **현재 분석 join 오류** | `HaidRaw_paIdx` ↔ `TracRaw_paIdx` | VECP와 PA라는 다른 index 공간을 숫자로 join |
| **추가 조사 필요** | `TracRaw_chi2NoVD`, 드문 invalid track 값 | native에서도 발견한 upstream 해석 문제가 남음 |
| **Native 오류가 재현되지 않음** | `Stic_*` 주요 운동량/에너지, Track perigee/matrix, beam spot | 별도 raw bank reader에서 올바른 위치/타입 사용 |

Producer 소스 기준은 [delphi-nanoaod `706d1dc…`의 raw writer](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp)와
[dstana `v270826`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car)이다.
확인에 사용한 cached raw writer가 이 revision의 파일과 byte 단위로 같음을 확인했다.
각 sample 생성 실행 바이너리의 build hash까지 추적한 것은 아니다.

## 2. Electron / muon ID: header를 물리량으로 저장

### Electron

| Branch | 전체 데이터에서 확인한 값 |
| --- | --- |
| `ElidRaw_tag` | 876,001개 전부 **25**: ELID module ID |
| `ElidRaw_gammaConvTag` | 876,001개 전부 **9040**: 실제 conversion tag가 아님 |
| `ElidRaw_refitMomentum.fCoordinates.fX` | 321,000–321,370: momentum이 아니라 packed descriptor |
| 대응 Y / Z | 각각 HPC 관련 값 / conversion radius를 읽음 |

Raw writer는 `LELID+1/+2/+3..5`를 읽는다.
Fortran `PSHELE`는 electron/conversion tag를 PA의 packed identification word에서 읽고,
refit momentum은 MAIN/ELTR에서 가져온다. ELID의 같은 위치를 tag와 momentum으로 해석하면 안 된다.
[Raw 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L588),
[Fortran `PSHELE`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

예: `ZKK/job_0`, raw ROOT entry 1, event 1의 두 electron-ID rows는
tag=25, conversion tag=9040, refit momentum=(321000, 0, 0)이다.
이 세 필드는 현재 값으로 electron PID나 conversion 선택에 사용할 수 없다.
`ElidRaw_paIdx`와 bank 존재 여부까지 잘못됐다는 뜻은 아니다.

### Muon

`MuidRaw_tag` 224,440개도 모두 **24**, 즉 MUID module ID다.
`MuidRaw_looseChi2`는 9–127인데 χ²가 아니라 block-presence pattern을 읽는다.
`MuidRaw_hitPattern`도 올바른 hit-pattern 위치가 아니다.

Fortran `PSHMUO`의 tag는 PA의 packed identification word에서,
very-loose χ²와 hit pattern은 MUID slot +4/+5에서 해당 block 존재 여부를 확인해 읽는다.
Raw writer는 +1/+2/+3을 사용한다.
[Raw 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L576),
[Fortran `PSHMUO`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

예: `Zmumu/job_0`, raw entry 2, event 2에서 tag=24, “looseChi2”=127이다.
세 물리량은 다시 저장해야 하며, `MuidRaw_paIdx`는 별개로 취급한다.

## 3. TPC packed 값: float의 비트를 integer로 읽음

대상은 `MtpcRaw_packedPadsSectors`, `MtpcRaw_packedWiresHits`다.
각 1,249,002개 값 모두 수십억 수준의 integer로 저장돼 있지만,
그 비트를 float32로 해석하면 정상적인 정수 packed 숫자가 나온다.

| 예시 | 저장 integer | 같은 비트의 float32 값 |
| --- | ---: | ---: |
| Pads/sectors | 1,195,722,496 | 50,503 |
| Wires/hits | 1,232,348,800 | 1,000,040 |

Producer가 `IQ(LMTPC+7/+8)`을 사용하지만 payload는 REAL에 숫자로 담긴 packed code다.
수치 해석에는 `Q`를 읽고 정수로 반올림해야 한다.
현재 저장 integer에 바로 `%100` 등의 packed-number 연산을 적용하면 잘못된 결과가 된다.
[Raw MTPC 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L1126),
[Fortran `GETDEDX`의 REAL packed 값 사용](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/dedxid.car).

이 두 필드의 문제를 모든 `MtpcRaw_*` 수치에 적용하지는 않았다.

## 4. FEMC shower의 layer metadata

`EmShower_detector=26`인 FEMC shower는 386,167개다.
그중 `EmShower_nLayers`가 0이 아닌 경우가 386,012개,
`EmShower_layerPattern`이 음수인 경우가 7,402개다.

FEMC의 기본 shower record는 energy와 position의 **4 words**다.
그런데 raw writer는 HPC와 동일하게 +5/+6을 읽어 nLayers/pattern을 저장한다.
이 값은 다른 record/footer를 읽은 것이므로 FEMC layer 정보가 아니다.
FEMC의 block 정보는 뒤의 별도 구조에서 읽어야 한다.
[Raw shower 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L338),
[Fortran `PSHEMC`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

예: `Zee/job_0`, raw entry 5, shower row 1:
energy=0.923738 GeV, nLayers=40, layerPattern=-36.

반면 HPC shower 887,792개에서는 layer 수/pattern 범위, bit count,
실제 emitted layer row 수가 모두 일치했다.
HPC layer 에너지 합이 shower 에너지와 1% 이상 다른 경우는 813,813개지만,
현재 reader 위치는 Fortran과 같다. 보정/packing 정의를 더 확인하기 전에는 이 차이를 저장 오류로 단정하지 않는다.

## 5. RICH: raw도 독립적인 정상 대조군이 아님

Raw RICH가 bank를 완전히 독립적으로 읽는 것은 아니다.
`HaidRaw_*` 저장도 native와 같은 SKELANA convenience arrays를 사용한다.

| Raw branch | 확인된 문제 |
| --- | --- |
| `HaidRaw_thetaGas` | gas angle이 아니라 `SIGG` resolution 영역을 읽음 |
| `HaidRaw_sigmaGas` | integer photon count의 float 비트; 0이 아닌 1,700,905개 모두 subnormal |
| `HaidRaw_nphGas` | observed가 아니라 expected count 영역 |
| `HaidRaw_nepGas` | expected count가 아니라 sector flag 영역 |
| `HaidRaw_flagGas` | gas COMMON 영역 밖 |
| `HaidRaw_nepLiq` | INTEGER를 float로 읽음; 0이 아닌 2,103,405개 모두 subnormal |
| `HaidRaw_flagLiq` | REAL을 integer로 읽음; 이번 데이터는 전부 0 |

Gas bank 길이 불일치와 liquid 타입 불일치는 native 보고서와 같은 원인이다.
Liquid angle/sigma/nph는 위치와 타입은 맞지만 이전 이벤트 값이 남는 문제가 있다.
[Raw RICH 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L677),
[Fortran COMMON 타입/길이](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/stdcdes.car).

예: `Zee/job_0`, raw entry 16, HaidRaw row 0의 liquid triplet은
angle=0.664330125, sigma=0.017584677, nph=16으로 직전 raw entry의 같은 row와 정확히 같다.
현재 raw quality와 다섯 dE/dx/RICH PID tags도 모두 0이다.
대응 native 이벤트에서도 현재 Haid 정보가 없는 charged Part에 같은 잔류값이 확인됐다.
모든 양수 RICH 값이 잔류값이라는 뜻은 아니다.
이벤트 초기화가 convenience arrays를 지우지 않는 원인은 raw에서도 그대로 적용된다.

### Raw에서 추가로 잘리는 quality word

`HaidRaw_richQuality`는 int8로 저장된다. Native의 대응 status word는 int32다.
전체 대응 rows 중 native 값이 int8 범위를 넘는 경우는 267,096개다.

- `Zee/job_0`, raw entry 3, row 0: native 320 → raw 64.
- `Zmumu/job_0`, raw entry 2, row 1: native 256 → raw 0.

따라서 raw quality=0을 “RICH 정보 없음”으로 판단하면 정보 손실을 놓칠 수 있다.
[int8 cast](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L675).

Native/raw quality가 다른 rows는 총 300,607개이고,
그중 39,759개는 native 값의 int8 cast와도 다르다.
각 writer는 별도 실행이므로 이 추가 차이의 원인을 truncation 하나로 설명하지 않았다.
RICH auxiliary arrays도 native와 항상 같지는 않으며, 처리 상태/설정 차이는 아직 확인되지 않았다.

## 6. 현재 raw RICH 분석의 index join 문제

**`HaidRaw_paIdx`라는 이름에도 불구하고 값은 PA index가 아니라 1-based VECP index다.**
반면 `TracRaw_paIdx`는 0-based PV/PA chain traversal index다.
Native `Part`도 VECP 순서를 사용하며, Fortran은 charged-first 재배열, PA 생략,
energy-recovery Part 추가 등을 수행할 수 있다.
Raw가 호출하는 `PSINI`의 기본값 `IFLODR=1`도 charged-first 순서를 설정한다.
따라서 `HaidRaw_paIdx - 1 == TracRaw_paIdx`는 오브젝트 동일성을 뜻하지 않는다.
[Raw Haid loop](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L662),
[Fortran `PSHCTK`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

전체 raw Track 1,438,959개 중 숫자 PA index가 native Part 범위 안인 1,438,929개를
같은 숫자의 Part row와 비교하면 charge가 다른 경우가 **357,629개, 24.85%**다.
나머지 30개는 native Part 수를 넘는다. 이는 raw index 범위 오류가 아니라 서로 다른 index 공간의 증거다.
Charge refit 등도 비교 차이에 기여할 수 있으므로 357,629개 전부를 단순 재배열 오류의 개수로 해석하지 않는다.

예: `Zee/job_0`, raw entry 1, event 1의 PA index 0은 raw charge=-1인데,
native Part row 0의 charge는 +1이다.

더 직접적인 예: `ZKK/job_0`, raw entry/event 390에는 raw charged Track 하나만 있고 PA ID는 1이다.
Charged-first VECP 순서에서는 해당 charged Haid row가 0인데,
현재 숫자 join은 Haid row 1을 선택한다. Liquid angle도 각각 0.6768022 / 0.67289835로 다르다.

`HaidRaw`는 2,291,622 rows로 native Part 수와 같고, writer는 모든 VECP rows를 저장한다.
Branch 설명의 “charged tracks”와 달리 neutral placeholders도 포함한다.
Row 수가 같다는 사실만으로 모든 raw/native 보조 값의 동일성을 보장할 수는 없다.

현재 [raw RICH `track_fields`](../python/03_rich/data.py)는 위 숫자 join을 사용한다.
**따라서 현재 raw RICH Track association/coverage 결과도 올바른 대응으로 검증되지 않았다.**
Producer에서 실제 PA↔VECP mapping을 함께 저장하거나 RICH rows에 실제 PA ID를 저장해야 한다.
Angular matching으로 이 내부 index 문제를 우회해 truth association이라고 해석해서는 안 된다.

## 7. Native와 달리 저장 위치가 맞는 주요 필드

### STIC

Raw `Stic_*`는 SSTC/MAIN bank를 직접 읽어 native의 `LENSTC` 오류를 우회한다.
25,864개 `Stic_energyFromMain`은 모두 양수이고 NaN/Inf/subnormal이 없다.
모든 row에서 MAIN vector의 크기와 에너지가 상대 차이 1e-6 이내로 일치했다.
`Stic_directionFromMain`은 unit direction이 아니라 **MAIN momentum vector**다.
[Raw STIC reader](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L508).

Raw EM/STIC PA ID overlap은 없었고, STIC PA에 charged raw Track이 연결된 경우도 없었다.
다만 “저장 위치와 수치가 정상”이라는 점이 STIC의 완전한 물리 효율을 입증하지는 않는다.

### Track perigee / weight matrix

Raw는 perigee의 d0/z0/θ/φ/1/R을 올바른 TRAC 위치에서 읽고,
`TracRaw_weightMatrix`도 올바른 15개 weight를 읽는다.
Native의 한-word offset 오류는 재현되지 않았다.
[Raw TRAC reader](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L773).

### Beam spot / simulated PV

`Event_beamSpot*`, `Event_simPV*`는 매 이벤트 채워진다.
모든 물리 이벤트에서 두 errorFlag는 0이며 beam spot sigma도 양수다.
Setup에서는 errorFlag=-1이다. Native의 frozen `BeamSpot_*`와 다르다.
[Raw beam spot 저장](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L1154).

## 8. 추가 조사 / 정보 부재 / 단위 설명

### Track fit 값

Raw Track 1,438,959개 중 `TracRaw_chi2NoVD<0`은 704,694개이고,
**100,041개는 ndfNoVD>0인데도 χ²가 음수**다.
Native와 selection/오브젝트 수가 달라 같은 count가 나와야 하는 것은 아니다.
Raw의 VD χ²에는 음수가 없었다.
`TracRaw_chi2VDHits`에는 최대 8.50e16의 큰 tail도 있어 그 해석 역시 추가 확인이 필요하다.

MAIN slot 16/17/18을 읽는 방식은 SKELANA와 같다.
FullDST MAIN 문서에는 이 영역이 momentum covariance로 기술된 경우가 있어 format 해석이 의심되지만,
이번 SDST의 실제 packing까지 확정한 것은 아니다. 음수를 모두 실패-fit sentinel로 처리하면 안 된다.
[PXDST MAIN layout](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/pxdst34.car).

추가 수치 이상은 θ>π인 Track 2개, invR=0인 Track 30개,
`MtpcRaw_zFitChi2<0`인 row 3개다. 원인은 확정하지 않았다.

### 비어 있는 정보

- `nTrackElement`: 모든 450,099 records에서 0. 전체 `TrackElement_*`가 비어 있어 이 데이터로 subdetector track-element study를 할 수 없다. Bank 부재와 reader lookup 문제 중 어느 쪽인지는 미확정.
- `VdAssocHit_signalToNoise`, `VdUnassocHit_signalToNoise`: 각각 3,757,652 / 13,939,050개 모두 0. 값이 0이라는 사실만으로 producer 오류로 분류하지 않았다.
- `GenPart_productionTime`, `GenPart_properLifetime`: 6,293,026개 모두 0. 이번 SDST의 `PSHLUJ`는 이 두 값을 명시적으로 0으로 채운다.
- `Event_fillNumber`: MC에서 fill metadata가 없어 모두 0.

### Gen vertex 단위 설명

Raw writer는 `GenPart_vertex`를 mm라고 설명하지만 이번 SDST reader는
simulated PV의 cm 좌표를 `VP(1..3)`에 스케일 변환 없이 복사한다.
같은 bank를 읽는 `Event_simPV*`도 cm다.

예: `ZKK/job_0`, raw entry 1의 Gen rows 6–10은
vertex=(-0.101905495, 0.000149826, -0.397978097)로 `Event_simPVXYZ`와 정확히 같다.
이번 SDST 값에 명목상 LUJETS mm convention을 적용해 10으로 나누면 잘못된 위치가 된다.
[Writer의 mm 설명](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp#L1266),
[Fortran `PSHLUJ`](https://gitlab.cern.ch/delphi/dstana/-/blob/v270826/src/car/skelana.car).

## 9. Index / truth 검사 결과와 현재 study 영향

전체 6,293,026 Gen rows에서 parent/firstChild/lastChild index의 범위와 sentinel은 정상이다.
Parent는 항상 child보다 앞에 있어 self-parent/cycle이 없다.
Parent가 있는 4,732,938 rows도 해당 parent's child range 안에 있었다.
대응 native/raw의 PDG/status/parent 배열은 모든 물리 이벤트에서 정확히 같았다.

아래 collection 내부 index는 전체 파일에서 유효한 0-based index였다.

| Branch | Rows |
| --- | ---: |
| `EmLayer_emShowerIdx` | 5,087,153 |
| `HadHit_hadShowerIdx` | 2,486,454 |
| `VdAssocHit_tracRawIdx` | 3,757,652 |
| `MtpcRaw_tracRawIdx` | 1,249,002 |

EM shower는 같은 PA에 여러 cluster가 있을 수 있다. 중복 PA 추가 rows 15,197개를
단순 중복 오브젝트로 삭제하면 안 된다.
EM shower 1,273,959개 중 722,179개는 charged raw Track의 PA에 연결된다.
따라서 raw EM shower 전체는 photon PID 후보와 같은 집합이 아니다.

Raw에는 `Part`/`SimPart`/`Photon`/`PhotonConv` collection 및 native Gen–Sim truth-link branches가 없다.
Gen ancestry 정의는 가능하지만 raw shower→Gen truth efficiency는 이 파일만으로 계산할 수 없다.
PA ID를 native Part index로 대입하는 것도 올바른 truth link가 아니다.

현재 분석에 대한 결론:

- **Raw Gamma:** energy/position 기반 EM/STIC angular study는 유지할 근거가 있다. FEMC layer metadata와 PID/내부 truth 해석에는 위 한계가 있다.
- **Raw Track:** perigee/matrix 저장 위치는 정상이나 zero curvature, 드문 invalid θ, NoVD χ²의 의미는 구분해야 한다.
- **Raw RICH:** branch payload와 Track join 문제가 모두 있어 현재 coverage/PID 성능으로 해석하기 어렵다.

우선 수정 대상은 raw ID bank reader, TPC packed 타입, FEMC metadata, RICH ABI/초기화/quality 타입과 PA↔VECP mapping이다.
이번 작업에서는 문제와 근거만 문서화했고 producer나 분석 결과를 수정하지 않았다.
