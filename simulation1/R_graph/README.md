# R_graph: 결과발표·논문용 그래프

결과발표회 원고의 그림을 만드는 R 코드입니다. 연구 목적 그림(그림 1~5)을 제외한 모든 그림은
시뮬레이션 결과 파일과 설정한 조건(`00_config.R`)으로 그립니다.

- 그래프에는 제목·설명을 넣지 않음 (캡션은 PPT·논문에서 직접 작성)
- 축·범례·패널 라벨은 한국어 (θ, P(θ), Bias, RMSE, GPCM·NRM 은 기호·약어 그대로)
- 글씨 기본 크기 20pt, 300 dpi PNG
- 비교되는 그래프는 축 범위와 패널 배열이 같음

## 실행

작업 디렉토리를 `simulation1/` 로 맞춘 뒤 실행합니다.

```r
setwd("경로/simulation1")
source("R_graph/run_all.R", encoding = "UTF-8")
```

- `01_collect_data.R` 이 결과 파일을 읽어 `<출력 폴더>/data/` 에 한 번만 저장합니다 (가장 오래 걸림).
  결과 파일이 바뀌지 않았으면 `run_all.R` 의 `COLLECT <- FALSE` 로 건너뛸 수 있습니다.
- 그림 하나만 다시 그리려면 해당 스크립트만 실행합니다. 예: `source("R_graph/06_fig_mechanism.R", encoding = "UTF-8")`
- 필요한 패키지: dplyr, tidyr, ggplot2, scales, ragg, car (없으면 자동 설치)

출력 폴더는 `00_config.R` 의 `NRM_SOURCE` 에 따라 나뉩니다 (서로 덮어쓰지 않음, `GRAPH_DIR` 로 직접 바꿀 수도 있음)
- `NRM_SOURCE <- "sim1"` (기본): `R_graph/output_sim1nrm/`
- `NRM_SOURCE <- "sim2"` (예전 방식): `R_graph/output/` (기존 그림이 있던 폴더)

각 출력 폴더 안
- 그림: `figures/*.png`
- 표: `tables/*.csv` (UTF-8 BOM, 엑셀에서 한글이 깨지지 않음)
- 정리 자료: `data/*.rds`

## 사용하는 자료

| 자료 | 경로 | 조건 |
|---|---|---|
| 시뮬레이션1 (GPCM 추정) | `output/estimated_params/`, `output/responses/` | 정규분포 27개 조건 (조건 코드 넷째 자리 2) |
| NRM 추정 (기본) | `output/nrm/estimated_params/` (`R/15_nrm_estimation.R` 출력) | 시뮬레이션1과 같은 응답, 정규분포 27개 조건 |
| NRM 추정 (예전 방식) | `../simulation2/output/nrm/estimated_params/` | `NRM_SOURCE <- "sim2"` 일 때만, 같은 27개 조건 |

NRM 자료는 `00_config.R` 의 `NRM_SOURCE` 로 고릅니다.

- `"sim1"` (기본): 시뮬레이션1 응답에 NRM 을 추정한 결과를 씁니다. 먼저 `source("R/15_nrm_estimation.R")` 로 추정해 두어야 합니다.
  GPCM 그림과 같은 응답·같은 조건이므로 경계모수 간격 1.5 도 그대로 씁니다.
- `"sim2"`: 예전 임시 방식. 시뮬레이션2에는 경계모수 간격 1.5 가 없어서 1.4 결과를 가져와 1.5 로 표시하고,
  편향·RMSE 의 참값은 1.4 로 계산합니다 (`SIM2_SUBSTITUTE_IV2`).

어느 쪽이든 정리 자료 이름은 `sim2_nrm_item1.rds` 이고, 그림 파일 이름(`fig12_sim2_*` 등)도 예전 그대로입니다.
NRM 결과가 없으면 NRM 그림만 건너뜁니다.

## 그림 목록

| 파일 | 원고 위치 | 내용 |
|---|---|---|
| `fig01_crc_discrimination` | 연구 목적 | 변별도 1 vs 2 의 범주반응곡선 |
| `fig02a_crc_boundary` | 연구 목적 | 범주반응곡선과 경계모수 위치 |
| `fig02b_boundary_curves` | 연구 목적 | 인접 두 범주만 남긴 경계반응곡선 |
| `fig03_crc_transposition` | 경계모수 전치 | 순서대로인 경계모수 vs 전치 (b3 > b4) |
| `fig04_scoring_function` | 경계모수 전치 | 채점함수 간격 1 / 0.66 / 0.33 (수직선) |
| `fig05_crc_scoring_function` | 경계모수 전치 | 채점함수 간격에 따른 범주반응곡선 |
| `fig05b_crc_interval` | 방법 (독립변인) | 경계모수 간격 0.5 / 1.0 / 1.5 에 따른 범주반응곡선 |
| `fig05c_crc_mean` | 방법 (독립변인) | 문항 심각도 0 / 0.5 / 1.0 에 따른 범주반응곡선 |
| `fig06_true_boundaries` | 방법 | 조건별 참 경계모수와 θ 분포 |
| `fig07_sim1_transposition` | 결과1 전치 비율 | GPCM 전치 비율 |
| `fig08_sim1_transposition_pairs` | 결과1 전치 비율 | GPCM 경계모수 쌍별 전치 비율 (누적 막대) |
| `fig09_sim1_bias_b3b4` / `_b1b4` | 결과1 Bias | GPCM Bias |
| `fig10_estimate_distribution` | 결과1 Bias | b̂3·b̂4 분포와 참값 |
| `fig11_1_b4_curves` | 결과1 편향 원리 | 범주 4·5: 참 / GPCM 문항특성곡선 (두 곡선 모두 참 b4 에 위치, 기울기만 다름) |
| `fig11_2_b4_curves_theta` | 결과1 편향 원리 | 11-1 + 범주 4·5 응답자 θ (y = 0, 반투명, 초록 θ < b4, 파랑 θ ≥ b4) |
| `fig11_3_b3_curves` | 결과1 편향 원리 | 범주 3·4: 참 / GPCM 문항특성곡선 (두 곡선 모두 참 b3 에 위치) |
| `fig11_4_b3_curves_theta` | 결과1 편향 원리 | 11-3 + 범주 3·4 응답자 θ |
| `extra_fig11_b3_above_*` | 질의응답 예비 | θ ≥ b3 응답자가 더 많은 조건의 11-3·11-4 (`FIG11_B3_ABOVE_COND`) |
| `fig12_sim2_transposition` | 결과2 전치 비율 | NRM 전치 비율 (그림 7과 같은 축) |
| `fig12b_transposition_difference` | 결과2 전치 비율 | NRM − GPCM 전치 비율 차이 |
| `fig13_sim2_rmse_b3b4` / `_b1b4` | 결과2 RMSE | NRM RMSE (이상치 제거 후) |
| `fig14_nrm_denominator` | 결과2 추가 설명 | 분모 â(ŝ5 − ŝ4) vs b̂4 |
| `fig15_nrm_cause` | 결과2 | NRM 전치 비율을 채점함수 전치 여부로 나눈 누적 막대 (그림 8과 같은 배치·축) |
| `fig16_summary_heatmap` | 논의 | 27개 조건 주요 결과 요약 |
| `extra_*` | 질의응답 예비 | GPCM RMSE, NRM Bias, NRM 쌍별 전치, KL 극한 vs 시뮬레이션 Bias, 채점함수 전치 사례 범주반응곡선 |

배치 규칙
- 전치 비율 막대 (그림 7·8·12·15): 9개 패널 — 행 = 채점함수 간격, 열 = 경계모수 간격, x = 문항 심각도
  y축은 0~100% 로 이 그림들과 예비 NRM 쌍별 그림이 모두 같음
- Bias·RMSE 선 그래프 (그림 9·13): 행 = 경계모수, 열 = 채점함수 간격, x = 문항 심각도, 색·모양·선 = 경계모수 간격
  Bias 그림끼리, RMSE 그림끼리 **모든 패널의 y축이 같음** (시뮬레이션1·2, b3b4·b1b4 판 모두)
- 그림 11-1~11-4 는 x·y 범위와 그림 크기가 모두 같아 슬라이드에서 겹쳐 넘겨도 곡선 위치가 같음

## 표 목록

| 파일 | 내용 |
|---|---|
| `tab_sim1_transposition_lrt` / `tab_sim2_transposition_lrt` | 전치 여부 로지스틱 회귀 우도비 검정 |
| `tab_sim1_bias_anova_b3` / `_b4` | GPCM (b̂ − b) Type III 분산분석 (η², 편 η²) |
| `tab_transposition_by_condition` | 조건별 전치 비율 (전체·쌍별) |
| `tab_error_by_condition` | 조건 × 경계모수별 Bias·RMSE |
| `tab_sim2_outliers` | NRM 조건별 이상치 제거 반복 수·비율 |
| `tab_nrm_cause_by_condition` | 조건별 NRM 전치 비율 (채점함수 전치 / 전치 없음) |
| `tab_fig11_summary` | 그림 11 에 쓴 값 (참값, 추정 평균, 대표 반복, 초록·파랑 점 수) |
| `tab_kl_limit` | KL 극한값과 시뮬레이션 평균 |

## 자주 바꿀 설정 (`00_config.R`)

| 설정 | 기본값 | 내용 |
|---|---|---|
| `SF_GAP_LABEL` | 1 / 0.66 / 0.33 | 채점함수 간격 라벨 (실제 4 − s4 = 1 / 0.67 / 0.34) |
| `ERR_CUT`, `OUTLIER_UNIT` | 3, `"rep"` | NRM 이상치: \|b̂ − b\| ≥ 3 인 값이 있는 반복 전체 제거 (simulation2 16번과 같음). `"param"` 이면 그 값만 제거 |
| `FIG11_B4_COND`, `FIG11_B3_COND` | 0.33 · 1.5 · 1.0 | 그림 11 조건 |
| `FIG11_B4_AUTO`, `FIG11_MIN_ABOVE` | TRUE, 3 | 그림 11-1·11-2 조건에 θ ≥ b4 응답자가 3명보다 적으면 b4 과소추정이 큰 다른 조건으로 자동 변경 |
| `FIG11_REP` | `NULL` | 그림 11 점에 쓸 반복 (NULL = 추정값이 조건 평균에 가장 가까운 반복) |
| `FIG11_POINT_ALPHA` | 0.25 | 그림 11 점 투명도 (겹칠수록 진해짐) |
| `TRANS_YMAX` | 1 | 전치 비율 그림의 y축 최댓값 (100%). 쌍별 누적 합이 넘으면 자동으로 늘림 |
| `FIG05_FIXED` | s4 = 3, 간격 1, 심각도 0 | 그림 5b·5c 에서 바꾸지 않는 변인의 값 |
| `FIG10_CONDS` | 2개 조건 | 그림 10 에 위아래로 놓을 조건 |
| `FIG14_COND`, `FIG14_YLIM` | 0.33 · 1.5 · 1.0, −10 ~ 15 | 그림 14 조건과 y 범위 (범위 밖 점 수는 콘솔에 출력) |
| `BASE_SIZE`, `DPI`, `SAVE_PDF` | 20, 300, FALSE | 글씨 크기, 해상도, PDF 동시 저장 |
| `FONT_FAMILY` | Windows: 맑은 고딕 | 한글 글꼴 |
| `LAB` | | 축·범례 라벨 |
