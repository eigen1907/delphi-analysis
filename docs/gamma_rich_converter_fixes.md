# Gamma / RICH converter fixes

2026-10-14, branch `fix/gamma-rich-converter`.
DELSIM/DELANA 코드는 변경하지 않았다. 기존 Z→KK PHOTOS-FSR 1,000개 이벤트의
SDST와 full FADANA를 수정된 NanoAOD writer로 다시 변환했다.

2026-10-15 코드 정리 후 두 converter를 다시 빌드하고, SDST/full DST 각각
100개 이벤트를 재변환했다. 기존 수정본과 모든 데이터 컬럼이 정확히 일치했으며,
native/raw RICH 11개 필드도 같은 PA 기준으로 일치했다.
Converter 커밋: Gamma `ff6051a`, RICH `47bc812`.

새 분석 입력:
`data/20261014_stic_validation/20261014_1kTest_ZKK_photosFSR_gamma_rich_fixed/final_root/`.
기존 Florian ROOT 파일은 수정 이전 producer의 결과이므로 아래 수정이 소급 적용되지 않는다.

## Gamma

| 수정 | 내용 |
| --- | --- |
| STIC | 이전 수정 유지: COMMON stride 9, silicon field 9. 실제 STIC 후보 25개 복구 |
| Raw FEMC | 4-word shower 밖의 layer metadata를 읽지 않음. 없는 `nLayers/layerPattern`은 `-1` |
| SimVtx | Mass code는 `KVTX(5)`, 좌표는 `QVTX(6..8)`에서 읽음 |
| Gen–Sim index | SH 순서의 `ISTSH/ISHST` 대신 LU/Gen 순서의 `ISTLU/ILUST` 사용 |
| PhotonConv | 저장된 reciprocal conversion vertex → neutral gamma parent Part → Sim gamma 연결로 `simPhotonIdx` 채움 |
| Raw electron/conversion | ELID header 대신 PA identification bits로 tag를 읽고, decoded electron bank에서 momentum 저장 |
| 설명 | `Photon_emEnergy`는 FEMC 전용; 현재 conversion type은 양수 21–24; DST Gen 위치는 cm |

`PhotonConv_simPhotonIdx`는 같은 vertex의 유효한 daughter와 reciprocal parent,
parent의 직접 Sim link 및 Sim gamma code 21이 확인되는 경우에만 채운다.
없는 연결은 `-1`로 남긴다. Angular matching을 사용하지 않는다.
기존 분석은 이미 parent/daughter ancestry를 사용하므로 matching 정의를 바꾸지 않았다.

검증 결과:

- Photon 1,280개, conversion 68개와 기존 운동량/Gen 분류 유지.
- Conversion의 직접 Sim-gamma index: 25개 유효, 43개 unresolved.
- 기존 Gen→Reco 개수, energy ratio, 유효성 및 detector flags는 1,000개 이벤트 모두 동일.
- SH/LU index 수정은 현재 표본의 저장값을 바꾸지 않음. 순서가 달라지는 입력에 대한 오류 예방.
- SimVtx 1,643개 좌표가 finite이며, Gen–Sim anchor의 위치 5,282쌍이 정확히 일치.
- FEMC 316개 layer metadata는 모두 `-1`; HPC 2,469개와 layer 12,775개의 기존 값 유지.
- Raw ELID 1,908개 tag/conversion tag/momentum이 같은 PA의 native 값과 정확히 일치.
  ELID 외 raw 컬럼 177개는 수정 전후 동일하여 기존 플롯 입력도 유지.

## RICH

| 수정 | 내용 |
| --- | --- |
| Gas COMMON | Bank 길이 6→5; expected-count auxiliary는 INTEGER, flag auxiliary는 REAL |
| Liquid COMMON | Expected-count auxiliary는 INTEGER, flag auxiliary는 REAL |
| 이벤트 처리 | `PSBEG`가 초기화한 decoded bank를 읽음. 별도 COMMON 초기화 불필요 |
| 값 저장 | 현재 이벤트의 decoded gas/liquid bank를 읽음. Fractional expected count 보존 |
| Raw quality | `int8` 대신 `int32`로 전체 quality word 저장 |
| Raw Track 연결 | `HaidRaw_paIdx`를 실제 0-based PA index로 저장 |
| Native/raw 대응 | 새 `Part_paIdx`와 `HaidRaw_paIdx`로 같은 PA를 대응. 없는 PA는 `-1` |

PA index 0은 유효하다. PA index와 VECP/Part row 번호는 같은 index 공간이 아니다.
1,000개 SDST 이벤트의 Part 5,388개 중 2,081개는 두 번호가 달랐다.

검증 결과:

| 검사 | SDST | Full FADANA |
| --- | ---: | ---: |
| 물리 이벤트 | 1,000 | 1,000 |
| PA로 대응한 native/raw RICH row | 5,388 | 5,090 |
| RICH 11개 필드 정확히 일치 | 100% | 100% |
| Subnormal payload | 0 | 0 |
| 잘못된 길이 / 중복 PA / Track 연결 누락 | 0 | 0 |

SDST native gas sigma의 subnormal 3,431개와 liquid expected-count의 subnormal
4,833개가 모두 사라졌다. 이전에 양수 angle이던 2,510개 row는 새 gas/liquid
bank가 모두 0이었다. 저장 오류/잔류 제거에 따른 변화이며 검출기 반응의 감소로 해석하지 않는다.

RICH 플롯의 coverage 정의는 **양수 normal Cherenkov angle이 저장된 Track의 비율**이다.
검증하지 않은 quality-bit cut이나 PID decision은 추가하지 않았다.
수정 이전 schema에서는 잘못된 값으로 coverage를 계산하지 않고 `00_data_validation`만 실행한다.

## 남은 제한

- STIC 25개에는 양방향 Sim link가 없다. Converter 수정으로 새 truth 출처를 만들지 않는다.
- `Photon`은 neutral EM-cluster 후보이며 photon PID를 보장하지 않는다.
- Conversion의 parent link가 없는 43개 직접 Sim index는 여전히 unresolved다.
- Gen time/lifetime은 해당 DST unpacker에서 제공되지 않는다.
- Raw EM shower reader는 SDST `EMNC` 형식용이다. Full FADANA의 빈 EmShower를 실제 shower 부재로 해석하지 않는다.
- Full FADANA의 raw ELID collection도 비어 있어 electron-ID 부재로 해석하지 않는다.
- Electron refit momentum에는 최대 약 714 GeV의 tail이 있다. 같은 값이 native에 저장되어 있으며,
  저장 위치 오류로 단정하지 않는다. 해당 fit 값의 물리 사용에는 별도 품질 검증이 필요하다.
- Full DST와 SDST의 RICH 값은 packing/반올림 차이가 있다. 각 입력의 native/raw 일치를 별도로 검증했다.

## 실행

```bash
uv run python scripts/02_gamma/plot.py -i data/20261014_stic_validation/20261014_1kTest_ZKK_photosFSR_gamma_rich_fixed
uv run python scripts/03_rich/plot.py -i data/20261014_stic_validation/20261014_1kTest_ZKK_photosFSR_gamma_rich_fixed
```

원본 입력과 재변환 파일은 `data/20261014_stic_validation/production/`에 보존했다.
