# 문서

실행 명령과 저장 구조는 [프로젝트 README](../README.md)를 따른다.

| 문서 | 내용 |
| --- | --- |
| [Branch 사용 규칙](branch_usage_summary.md) | Native/raw의 역할, index와 결측값 해석 |
| [Native branch 사전](nanoaod_branch_audit.md) | `nanoaod.root`의 확인된 문제와 사용 범위 |
| [Raw SDST branch 사전](nanoaod_raw_sdst_branch_audit.md) | `nanoaod_raw_sdst.root`의 확인된 문제와 사용 범위 |
| [Converter 수정](gamma_rich_converter_fixes.md) | Gamma/RICH 디버깅으로 바뀐 ROOT branch |
| [Photon fiducial 정의](fiducial_regions.md) | Gen/Reco 선택, detector 영역, 효율의 분자·분모 |
| [BDT 입력](bdt_input_branches.md) | 현재 feature 4종과 공통 입력 |

원본 ROOT의 tree는 `Events`다. 준비된 dataset은 같은 이벤트를 정렬해
`analysisLevel/Events`와 `raw_sdst/Events`에 함께 저장한다.

`20260828_florian`은 수정 전 converter의 출력이다.
`20261015_jongwon`에는 [수정된 Gamma/RICH schema](gamma_rich_converter_fixes.md)가 적용됐다.
기존 Florian 파일의 값이 코드 수정만으로 바뀌지는 않는다.

전체 수치 감사 자료는 [native CSV](nanoaod_numeric_scan.csv)와
[raw CSV](nanoaod_raw_sdst_numeric_scan.csv)에 보존한다. 두 CSV는 Florian 원본 기준이다.
