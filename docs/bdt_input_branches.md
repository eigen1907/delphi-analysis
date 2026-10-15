# BDT 입력

분류 대상은 `Zee=0`, `Zmumu=1`, `ZKK=2`, `Zpipi=3`, `Ztautau=4`다.
`train`으로 학습하고 `val`로 early stopping, `test`로 평가한다.

## Feature 설정

| 설정 | RAW 입력 | 최종 feature 수 |
| --- | --- | ---: |
| [detector](../config/pid/features/detector.json) | 아래 detector collection | 743 |
| [detector_pid](../config/pid/features/detector_pid.json) | detector + `MuidRaw`, `ElidRaw`, `HaidRaw` PID | 856 |
| [detector_rich](../config/pid/features/detector_rich.json) | detector + RICH 8개 | 823 |
| [detector_pid_rich](../config/pid/features/detector_pid_rich.json) | detector_pid + RICH 8개 | 936 |

기본은 `detector_pid`. 모든 설정에 아래 analysis-level 5개가 포함된다.
전체 branch 선택은 각 JSON, 최종 열 순서는 저장 모델의 `features`에 기록된다.

`raw_sdst/Events`의 detector collection:
`Event`, `TracRaw`, `EmShower`, `EmLayer`, `HadShower`, `HadHit`,
`Stic`, `MtpcRaw`, `VdAssocHit`, `VdUnassocHit`, `Vtx`.

RICH 8개:
`HaidRaw_thetaGas`, `HaidRaw_sigmaGas`, `HaidRaw_nphGas`, `HaidRaw_nepGas`,
`HaidRaw_thetaLiq`, `HaidRaw_sigmaLiq`, `HaidRaw_nphLiq`, `HaidRaw_nepLiq`.
`detector_rich`에는 이 8개만 추가한다. `detector_pid_rich`는 기존 `nHaidRaw`와 PID tags도 유지한다.
Quality/sector flag는 추가하지 않는다.

`analysisLevel/Events`의 공통 5개:
`Event_chargedMult`, `Event_neutralMult`, `Event_totalChargedEnergy`,
`Event_totalEMEnergy`, `Event_totalHadronicEnergy`.
Gen/Sim truth, event/job ID, label은 모델 입력에 넣지 않는다.

## 추출·학습 규칙

- Track은 `abs(TracRaw_invR)` 오름차순의 상위 10개를 사용한다. 공분산은 상위 4개 × 15개 성분이다.
- Track 관련 collection은 저장된 PA/track index로 연결한다. 두 tree의 입자 row를 같은 번호로 연결하지 않는다.
- 부족한 슬롯은 NaN이다. Florian의 연결 불가능한 HaID vector도 NaN이며 선택된 collection count는 유지한다. 따라서 Florian의 추가 RICH는 새 측정 정보를 제공하지 않는다.
- 학습 설정은 [default.json](../config/pid/hyperparameters/default.json) 하나다: 최대 70 boosting rounds, tree당 최대 5 leaves. 기본 seed 1907, CPU 32.

구현: [feature 추출](../src/event_classification/bdt.py), [dataset 로더](../src/event_classification/dataset.py).
실행 명령은 [README](../README.md)를 따른다.
