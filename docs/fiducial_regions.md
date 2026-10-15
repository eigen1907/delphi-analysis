# Photon fiducial 정의

기준 구현은 [data.py](../src/sample_validation/gamma/data.py),
[angular matching](../src/sample_validation/gamma/angular_matching.py),
[detector plots](../src/sample_validation/gamma/detector_plots.py)다.
기본 입력은 `output/dataset/20261015_jongwon`의 **train**, 기본 stage는 `all`이다.

## Gen 선택

`GenPart_pdgId=22`, `GenPart_status=1`인 stable photon에서 시작한다.
Gen ancestry로 beam ISR, non-beam ISR, FSR, 나머지(`decayed`)를 분류한다.
Beam ISR은 ISR 중 `px=py=0`이고 parent가 incoming electron인 경우다.

| Matching population | 선택 |
| --- | --- |
| `gen_gamma` | Beam ISR을 제외하고 아래 EM fiducial을 통과 |
| `gen_isr` | Non-beam ISR이며 EM fiducial을 통과 |
| `gen_no_isr` | FSR 또는 `decayed`이며 EM fiducial을 통과 |

`01_sample_distribution`의 Gen 분포는 beam ISR과 영역 밖 photon도 포함한다.
아래 matching 효율의 분모와 구분해야 한다.

## EM 영역

θ는 **Gen photon 방향의 polar angle**, 단위는 degree다. 모든 경계는 엄격한 부등호다.

| 영역 | θ 범위 |
| --- | --- |
| HPC | `44 < θ < 88` 또는 `92 < θ < 136` |
| FEMC | `11 < θ < 36` 또는 `144 < θ < 169` |
| STIC | `2 < θ < 10` 또는 `170 < θ < 178` |

EM fiducial은 세 영역의 합집합이다. 경계와 사이 gap은 제외한다.
그림의 **All ranges (0–180°)**도 이미 fiducial을 통과한 Gen의 합이며, 4π 분모가 아니다.
Reco 후보에 같은 θ cut을 미리 적용하지 않는다.

## Reco와 matching

| 소스 | Reco 정의 |
| --- | --- |
| Native `Photon` | 저장된 neutral EM-cluster 후보 전체. Photon PID 확정 집합이 아님 |
| Native `Photon_PhotonConv` | Photon + conversion. 저장된 parent/daughter 관계로 중복 표현을 정리 |
| Native `Part_pdgId_21` | DELPHI mass code 21인 Part. 미식별 neutral의 기본 code도 포함 가능 |
| Raw `EmShower`, `Stic`, `EmShower_Stic` | 유한한 양수 에너지와 유한한 nonzero 방향을 가진 cluster. Charged PA cluster도 유지 |

Raw EM 방향은 `cluster position − simulated PV`, STIC은 저장 MAIN momentum으로 만든다.
Raw에는 Reco→Gen Sim link가 없어 angular study만 한다.

- **Angular:** Gen마다 독립적인 3D opening-angle cone을 사용한다. 한 Reco가 여러 Gen cone에 들어갈 수 있다. One-to-one matching이 아니다.
- **Truth link:** Native의 양방향 저장 관계를 비교해 유일한 Gen origin을 사용한다. Conflict/no-origin을 angular로 대신 채우지 않는다. Sim ancestry는 첫 Gen anchor에서 멈춘다.
- Angular cut은 `0.25, 0.5, 1, 2, 3, 5, 10, 20°`를 scan한다. Native collection별로 전체 샘플의 non-beam ISR truth 성공 중 99%를 유지하는 최소 cut을 고른다. 이 기준은 fiducial 적용 전이다. 만족하는 cut이 없으면 20°다.
- Raw는 Native `Photon_PhotonConv`에서 얻은 cut을 사용한다. 코드의 초기 참고값 3°를 고정 default cut으로 해석하지 않는다.

## Regional과 flag의 분자·분모

에너지 구간·population·Gen 영역을 통과한 photon 집합을 `G`라 하면:

| 그림 | 분모 | 분자 |
| --- | --- | --- |
| `eff_det_regional` | `G`의 모든 Gen photon | Matching된 Reco가 하나 이상 있는 Gen |
| `eff_det_flag` | 동일한 `G` | Matching된 Reco 중 해당 detector tag가 있는 Gen |

All flag는 HPC/FEMC/STIC 중 하나라도 있으면 성공이다.
각 detector bin은 Gen θ로 나누며, regional 분자는 Reco tag를 요구하지 않는다.
미검출 Gen도 분모에 남기고, 분모가 0인 점은 표시하지 않는다.

Native는 저장된 `Photon_isHPC/isFEMC/isSTIC`를 그대로 사용한다.
Conversion에는 detector tag를 만들지 않는다. Raw tag는 EM detector ID 9=HPC, 26=FEMC 및 STIC collection으로 정한다.
따라서 conversion만 복구한 Gen은 regional 성공이어도 flag 성공이 아닐 수 있다.
Raw 결과는 photon PID 효율 대신 **angular EM coverage**로 표시한다.

## 에너지 구간과 출력

에너지는 **Gen photon의 E [GeV]**다. 같은 구간을 분자·분모에 적용하며 Reco pool은 바꾸지 않는다.

| 구간 | Gen E | 파일 suffix |
| --- | --- | --- |
| All, 기본 | `0 ≤ E`, 상한 없음 | 없음 |
| Low | `0 ≤ E < 2` | `_energy0-2` |
| Mid | `2 ≤ E < 10` | `_energy2-10` |
| High | `10 ≤ E` | `_energy10-` |

입력 Gen four-vector는 유한하고 momentum과 E가 양수여야 한다.
`eff_det_flag.png`, `eff_det_regional.png`와 세 구간별 파일을 함께 만든다.
같은 구간은 `eff_vs_cos_theta`와 `event_energy_ratio`에도 적용된다.

기본 출력:
`output/plots/<campaign>/03_sample_validation/train/02_gamma/<stage>/matching_result/<sample>/<source>/<collection>/`.
에너지 합은 unmatched=0을 포함한다. 겹친 Reco 표현의 에너지가 모호하면 NaN으로 남기며,
angular event 합은 여러 Gen cone에 들어간 Reco를 한 번만 센다.

실행: `bash scripts/03_sample_validation/02_gamma/plot.sh`.
`--input`, `--split`, `--stage`로 대상을 바꿀 수 있다.
