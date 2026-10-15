# DELPHI Analysis

DELPHI 샘플 생성, ROOT 검증, BDT 이벤트 분류 코드.

## 설치와 구조

Python 3.12 환경에서 실행한다.

```bash
uv sync --locked
```

```text
scripts/
  01_sample_generation/       생성·시뮬레이션·변환
  02_sample_preparation/      train/val/test 준비
  03_sample_validation/       00_overall, 01_track, 02_gamma, 03_rich
  04_event_classification/    학습·평가·추론
src/                         재사용하는 Python 코드
config/pid/                  BDT feature·hyperparameter JSON
docs/                       branch와 분석 정의
data/                        원본 샘플
output/
  dataset/                   paired ROOT dataset
  plots/                     검증·BDT 그림
  models/                    모델·평가 지표
```

`scripts/`의 Bash recipe에 기본 옵션을 적었다. 옵션을 수정하거나 명령 뒤에
인자를 붙여 덮어쓸 수 있다. 상대 경로는 저장소 루트 기준이다.
Python 패키지는 `src/`에서 설치되므로 `PYTHONPATH` 설정은 필요 없다.
데이터·출력·컨테이너·캐시는 Git에서 제외한다.

## 1. 샘플 생성

Apptainer와 준비된 SIF 이미지가 필요하다.
[설치·이미지 준비·중간 결과](scripts/01_sample_generation/README.md)를 참고한다.

```bash
mkdir -p data/20261015_jongwon
nohup bash scripts/01_sample_generation/run.sh --sample all --ncpu 100 \
    > data/20261015_jongwon/run.log 2>&1 &
```

5개 채널 각각 5,000 이벤트 × 20 jobs를 실행한다. `--ncpu`는 전체 동시 job 수다.
DELSIM/DELANA 결과와 단계별 로그는 각 샘플의 `production/job_N/`에 보존한다.

## 2. Dataset 준비

각 원본 job의 `nanoaod.root`와 `nanoaod_raw_sdst.root`를 run/event로 맞추고,
유효한 PHOTOS 이벤트를 중복 없이 train/val/test로 나눈다. 원본은 유지한다.

```bash
bash scripts/02_sample_preparation/prepare_dataset.sh \
    --input data/20261015_jongwon --output output/dataset/20261015_jongwon \
    --counts 60000 15000 15000 --seed 1907

bash scripts/02_sample_preparation/prepare_dataset.sh \
    --input data/20260828_florian --output output/dataset/20260828_florian \
    --samples Zee ZKK Zpipi Ztautau --counts 60000 15000 15000 --seed 1907
bash scripts/02_sample_preparation/prepare_dataset.sh \
    --input data/20260828_florian --output output/dataset/20260828_florian \
    --samples Zmumu --counts 60000 15000 14999 --seed 1907
```

`--counts` 순서는 TRAIN VAL TEST다. Florian Zmumu는 원본 이벤트 1개가 없어
test만 14,999개다. 기존 dataset 파일은 덮어쓰지 않는다.

```text
output/dataset/<campaign>/
  metadata.json
  <sample>_photosFSR/{train,val,test}.root
    analysisLevel/Events
    raw_sdst/Events
```

두 tree는 `(job_id, run, event)` 순서가 같다. 기존 branch·중첩 구조를 유지하고
`job_id`를 추가한다. 분할은 이벤트 단위이며 split 사이에 중복은 없다.
Metadata에 원본 경로, 제외한 job/row와 실제 샘플별 이벤트 수를 기록한다.
별도 ROOT 설치 없이 uproot으로 처리한다.

## 3. Train dataset 검증

```bash
bash scripts/03_sample_validation/plot.sh --input output/dataset/20261015_jongwon
bash scripts/03_sample_validation/plot.sh --input output/dataset/20260828_florian
```

기본 split은 `train`이다. `--split val` 또는 `--split test`로 변경한다.
각 분석의 `plot.sh`는 전체 단계를, 번호가 붙은 recipe는 해당 단계만 실행한다.
`--output`을 생략하면 다음 경로를 사용한다.

```text
output/plots/<campaign>/03_sample_validation/<split>/
  00_overall/{00_branch_distribution,01_source_comparison,02_gen_distribution,
              03_gen_source_comparison,04_reco_distribution}/
  01_track/{00_data_validation,01_sample_distribution,02_angular_matching,03_truth_link_matching}/
  02_gamma/{00_data_validation,01_sample_distribution,02_angular_matching,03_truth_link_matching}/
  03_rich/{00_data_validation,01_sample_distribution,02_angular_matching,03_truth_link_matching}/
```

Photon 효율은 [fiducial region·matching·에너지 구간](docs/fiducial_regions.md)을 따른다.
이전 RICH schema는 유효한 트랙 대응을 보장하지 않아 저장 상태 검사만 수행한다.

## 4. BDT 학습·평가·추론

`CAMPAIGN`은 두 dataset 중 하나를, `FEATURE_SET`은 `detector`, `detector_pid`
(기본), `detector_rich`, `detector_pid_rich` 중 하나를 선택한다.
고정 hyperparameter는 `default`이며, 5종 `Zee`, `Zmumu`, `Ztautau`, `ZKK`,
`Zpipi`를 분류한다. [입력 feature](docs/bdt_input_branches.md)를 참고한다.

```bash
CAMPAIGN=20261015_jongwon FEATURE_SET=detector_pid bash scripts/04_event_classification/01_model_train.sh
CAMPAIGN=20261015_jongwon FEATURE_SET=detector_pid bash scripts/04_event_classification/02_model_evaluate.sh
CAMPAIGN=20261015_jongwon FEATURE_SET=detector_pid bash scripts/04_event_classification/03_model_apply.sh
```

Train으로 학습하고 val로 early stopping한다. Test는 평가에만 사용한다.

```text
output/models/<campaign>_<feature>_default/                    모델·metrics·importance
output/plots/<campaign>/04_event_classification/
  02_model_evaluate/<feature>_default/                        confusion matrix·importance
output/dataset/<campaign>/predictions/<feature>_default/      추론 ROOT
```

현재 8개 실험의 요약은 `output/models/results.csv`에 있다.
Branch 해석·변환기 수정 내역은 [문서 목차](docs/README.md)를 참고한다.
