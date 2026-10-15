# Raw SDST branch 사전

대상: `20260828_florian` PHOTOS-FSR 원본 `nanoaod_raw_sdst.root:Events`.
2026-10-14에 100파일의 450,099 records를 검사했다.
물리 이벤트는 native와 같은 449,999개이며, 파일마다 setup record 1개가 추가된다.
185개 필드 중 수치 payload 171개, collection leaf 길이 135개를 확인했다.
원자료는 [CSV](nanoaod_raw_sdst_numeric_scan.csv)에 보존한다.

현재 dataset에서는 `raw_sdst/Events`로 읽는다.
아래 “수정”은 [converter 변경](gamma_rich_converter_fixes.md)을 뜻하며 **기존 Florian 파일에는 소급되지 않는다.**

## 주요 branch

| Branch / collection | Florian에서 확인한 내용 | 현재 사용·수정 상태 |
| --- | --- | --- |
| `GenPart_*` | Native와 PDG/status/parent가 모든 물리 이벤트에서 일치. Ancestry index 범위 정상 | 위치는 cm. 시간/lifetime 없음 |
| `ElidRaw_tag`, `ElidRaw_gammaConvTag`, `ElidRaw_refitMomentum` | Tag=25, conversion tag=9040 등 ELID header/다른 값을 잘못 저장 | 수정 producer에서 PA ID bits와 decoded electron momentum으로 교정 |
| `MuidRaw_tag`, `MuidRaw_looseChi2`, `MuidRaw_hitPattern` | Tag=24(module ID), χ²는 block pattern, hit pattern도 잘못된 위치 | 이번 Gamma/RICH 수정 대상 아님. Bank 존재·PA index와 실제 PID 값의 해석을 구분 |
| `MtpcRaw_packedPadsSectors`, `MtpcRaw_packedWiresHits` | REAL packed 숫자를 integer 비트로 저장. 예: 1,195,722,496의 float 비트는 50,503 | 이번 수정 대상 아님. 저장 integer에 바로 `%100` 등을 적용하지 않음. 다른 dE/dx 측정값까지 같은 오류는 아님 |
| `EmShower_nLayers`, `EmShower_layerPattern` | FEMC 4-word record 밖의 +5/+6을 읽음 | 수정 producer는 FEMC에 `-1`. HPC는 기존 값 유지 |
| `EmLayer_*` | HPC layer 수/index/bit pattern 일치 | Layer 에너지 합과 보정 shower 에너지의 차이를 저장 오류로 단정하지 않음 |
| `HaidRaw_thetaGas/sigmaGas/nphGas/nepGas/flagGas` | Gas 길이 오류. `sigmaGas` nonzero 1,700,905개 모두 subnormal | 수정 producer에서 decoded gas bank로 교정 |
| `HaidRaw_thetaLiq/sigmaLiq/nphLiq/nepLiq/flagLiq` | `nepLiq` nonzero 2,103,405개 모두 subnormal; flag 타입 오류·이벤트 잔류 | 수정 producer에서 decoded liquid bank로 교정 |
| `HaidRaw_richQuality` | int8 축소. Native 256이 raw 0이 되는 등 정보 손실 | 수정 producer는 int32. 추가 writer 차이도 있어 과거 불일치를 truncation 하나로 설명하지 않음 |
| `HaidRaw_paIdx`, `nHaidRaw` | 1-based VECP 번호이며 neutral placeholder 포함 | 수정 producer는 실제 0-based PA ID/`-1`. Row 수가 같아도 Part/PA index를 같은 숫자로 연결하지 않음 |
| `Stic_*` | MAIN 에너지/momentum 저장 위치 정상. 25,864개 모두 양수·finite, momentum 크기와 에너지 일치 | `Stic_directionFromMain`은 momentum(GeV). 검출 효율 전체가 검증됐다는 뜻은 아님 |
| `TracRaw_*` | Perigee 5개·weight 15개 위치 정상. Native offset 오류 없음 | 드문 invalid θ/zero curvature와 fit-quality 의미는 구분 |
| `TracRaw_chi2NoVD`, `TracRaw_chi2VDHits`, `MtpcRaw_zFitChi2` | 양수 ndf에서도 음수 NoVD χ² 100,041개; 큰 VD-hit χ² tail, z-fit χ² 음수 3개 | 실패-fit sentinel로 일괄 처리하지 않음. Packing/상위 정의 추가 확인 필요 |
| `Event_beamSpot*`, `Event_simPV*` | 이벤트마다 갱신. 물리 이벤트 errorFlag=0, setup은 −1 | Native의 초기값 `BeamSpot_*`와 구분 |
| `HadShower`, `HadHit`, `Vtx` | 일부 저장 위치/index 확인 | 전체 물리량·fit 성능 검증 완료를 뜻하지 않음 |

## 빈 정보와 결측

| Branch | Florian에서 확인한 제한 |
| --- | --- |
| `TrackElement_*`, `nTrackElement` | 모든 record에서 빈 collection. Bank 부재/lookup 원인은 미확정 |
| `VdAssocHit_signalToNoise`, `VdUnassocHit_signalToNoise` | 전부 0. 기본 hit/index까지 무효라는 뜻은 아님 |
| `GenPart_productionTime`, `GenPart_properLifetime` | SDST unpacker가 0으로 채움 |
| `Event_fillNumber` | MC fill metadata 없음 |

Setup은 `nGenPart=0`, event number=0이다. 이 record 제외와 photon/track이 0개인
물리 이벤트를 버리는 것은 다르다. 드문 track θ>π 2개, invR=0 30개도 확인됐다.

## 연결과 photon 해석

- `EmLayer_emShowerIdx`, `HadHit_hadShowerIdx`, `VdAssocHit_tracRawIdx`, `MtpcRaw_tracRawIdx`는 검사한 전체 파일에서 유효한 0-based index였다.
- 한 PA에 여러 EM shower가 있을 수 있다. 같은 PA라는 이유로 cluster를 중복 삭제하지 않는다.
- EM shower 1,273,959개 중 722,179개가 charged raw track PA에 연결됐다. 전체 EM shower를 photon PID 후보로 해석하지 않는다.
- Raw에는 native `Part`/`SimPart`/`SimVtx`/`Photon`/`PhotonConv`와 Reco→Gen truth link가 없다. Gen ancestry와 shower angular coverage는 계산할 수 있다.
- 수정 전 PA와 Part 번호를 그대로 비교하면 charge 차이도 있었다. 이 숫자 join이나 angular matching을 내부 truth 연결로 대신 쓰지 않는다.
- RICH quality=0을 정보 부재로 단정하지 않는다. Florian은 [producer 재변환](gamma_rich_converter_fixes.md) 전까지 RICH response/coverage 해석을 보류한다.

근거: [수치 CSV](nanoaod_raw_sdst_numeric_scan.csv),
[당시 raw writer `706d1dc`](https://github.com/jingyucms/delphi-nanoaod/blob/706d1dcab7c5105e27087ff3aad1d6a8e8087e6a/delphi-raw-nanoaod/src/raw_nanoaod_writer.cpp),
[공통 사용 규칙](branch_usage_summary.md). 원본 실행 바이너리의 build hash까지 확인한 감사는 아니다.
