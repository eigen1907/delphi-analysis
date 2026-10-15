# Branch 사용 규칙

원본과 dataset의 대응은 다음과 같다.

| 원본 | Dataset tree | 주 용도 |
| --- | --- | --- |
| `nanoaod.root:Events` | `analysisLevel/Events` | 재구성 Part/Photon, Gen–Sim–Reco 관계 |
| `nanoaod_raw_sdst.root:Events` | `raw_sdst/Events` | Detector 측정값, raw track/PID, calorimeter cluster |

## 현재 사용할 때

| 용도 | 입력 | 해석·제한 |
| --- | --- | --- |
| Gen/ISR 분류 | 두 소스의 `GenPart_*` | PDG, status, ancestry 사용. 시간/lifetime 정보는 없음 |
| Native track 운동량 | `Trac_partIdx → Part_fourMomentum` | Native perigee/matrix의 offset 문제와 구분 |
| Raw track | `TracRaw_*` | Perigee/matrix 저장 위치 정상. NoVD χ²의 음수·드문 invalid track은 별도 확인 |
| Native photon | `Photon`, `PhotonConv`, `Part`, 저장된 truth 관계 | EM 후보와 실제 photon PID를 구분. [선택·효율 정의](fiducial_regions.md) 참고 |
| Raw photon study | `EmShower`, `Stic` | Charged cluster도 포함한 angular coverage. Raw만으로 Reco→Gen truth 효율 계산 불가 |
| Jongwon RICH | `HaidRaw_*` 또는 `Rich_*` | 수정된 PA 연결 사용. Quality bit의 의미·cut은 미검증 |
| Florian RICH | 기존 `HaidRaw_*`, `Rich_*` | ABI·잔류값·PA 연결 문제. 물리 response/coverage 해석 보류 |
| Vertex/부가 정보 | [Native 사전](nanoaod_branch_audit.md), [Raw 사전](nanoaod_raw_sdst_branch_audit.md) | 초기값·빈 payload·버전 의존 정보를 구분 |

“저장 위치 정상”은 검출기 성능이나 PID 성능의 검증 완료를 뜻하지 않는다.
현재 BDT의 raw χ² 등은 저장값을 feature로 사용하며, 물리적인 fit-quality cut의 검증과는 별개다.

## Index와 이벤트

- 두 소스는 `(job_id, run, event)`로 대응한다. 행 번호로 원본 파일을 붙이지 않는다.
- Raw의 빈 setup record는 물리 이벤트가 아니다. Dataset 생성에서 제외한다.
- PA index와 native Part/VECP row는 다른 공간이다. 수정된 `Part_paIdx`와 `HaidRaw_paIdx`는 0-based이고 없는 연결은 `-1`이다. **PA 0은 유효하다.**
- 수정 전 Florian의 `HaidRaw_paIdx`는 1-based VECP 번호다. 단순히 1을 빼서 `TracRaw_paIdx`에 연결할 수 없다.
- Native `Trac`/`Haid`/`Rich`는 neutral placeholder를 포함한 Part 순서를 따른다.
- Native `Vtx_incomingIdx/outgoingIdx`는 1-based/0 sentinel, `SimVtx_*` index는 0-based/`-1` sentinel이다.

## 값과 truth 관계

- `GenPart_pdgId`는 표준 PDG다. Native `Part_pdgId`, `SimPart_pdgId`는 DELPHI mass code다: gamma=21, e−=2, e+=−2.
- Forward/reverse truth가 한 Gen에만 연결되면 사용하고, conflict/no-origin을 억지로 채우지 않는다. Sim ancestry는 첫 유효 Gen anchor에서 멈춘다.
- All-zero, `-1`, 빈 collection은 branch별 정의로 해석한다. `isfinite`만으로 잘못 읽은 값을 걸러낼 수 없다.
- CSV의 `values`는 펼친 숫자 수다. `subnormal`은 0이 아닌 float의 절댓값이 최소 normal 값보다 작은 경우다.
- `Photon_emEnergy`는 FEMC 에너지다. HPC 에너지는 별도 branch에 있다. Raw `Stic_directionFromMain`은 단위벡터가 아닌 momentum(GeV)이다.
- `Event_totalChargedEnergy`는 선택된 charged energy, `Event_totalEMEnergy`는 neutral EM energy다. 전체 Part 에너지 합과 같지 않을 수 있다.

디버깅으로 실제 바뀐 값과 설명만 고친 필드는 [Converter 수정](gamma_rich_converter_fixes.md)에 구분했다.
