# 시뮬레이션 1: GPCM 경계 모수 전치 연구

## 연구 목적

GPCM의 채점함수(재모수화된 NRM)가 등간격이 아님을 무시하고 문항모수를 추정하는 절차가 경계 모수의 전치를 유발하는지 검증한다.

---

## 파일 구조

```
simulation1/
├── R/
│   ├── 00_setup.R           패키지, 전역 상수, rsn 보정값 계산
│   ├── 01_seed_generator.R  전체 시드 테이블 생성 (최초 1회만 실행)
│   ├── 02_item_params.R     조건 코드 파싱 + 문항 모수 생성
│   ├── 03_response.R        GPCM 확률 계산 + 가상 응답 생성
│   ├── 04_estimation.R      mirt GPCM 모수 추정
│   ├── 05_logger.R          로그 유틸리티
│   └── 06_main.R            메인 실행 스크립트 (조건별 단위 실행)
├── output/
│   ├── seed_table.rds            전체 시드 테이블 (108조건 × 1000반복)
│   ├── progress/
│   │   └── cond[C]_completed.txt 완료된 반복 번호 목록 (재개용)
│   ├── logs/
│   │   ├── temp/                 병렬 실행 시 임시 로그 (자동 삭제됨)
│   │   ├── gen_log_cond[C].txt   데이터 생성 로그
│   │   └── est_log_cond[C].txt   모수 추정 로그
│   ├── true_params/
│   │   └── YYYYMMDD_seed[S]_cond[C]_rep[R]_item_params.rds
│   ├── responses/
│   │   └── YYYYMMDD_seed[S]_cond[C]_rep[R]_response.rds
│   └── estimated_params/
│       └── YYYYMMDD_seed[S]_cond[C]_rep[R]_est_params.rds
└── README.md
```

파일명 형식: `YYYYMMDD_cond[C]_rep[R]_seed[S]_[type].csv`

파일명 변수 의미:
- `[C]`: 조건 코드 (예: `1233`)
- `[R]`: 반복 번호 (4자리 0패딩, 예: `0001`)
- `[S]`: 해당 반복의 시드 값

---

## 설계

### 독립변수 (108 조건 = 3 × 3 × 3 × 4)

| 자리 | 독립변수 | 수준1 | 수준2 | 수준3 | 수준4 |
|------|---------|-------|-------|-------|-------|
| 1 | 채점함수 sf[4] 값 | 3 | 3.33 | 3.66 | — |
| 2 | 경계모수 간격 | 2 | 1.5 | 1 | — |
| 3 | 문항 심각도(경계모수 평균) | 0 | 1.5 | 3 | — |
| 4 | 능력모수 분포 | 정적편포 | 정규분포 | 부적편포 | 균등분포 |

**조건 코드 예시:** `1233` → 채점함수 간격=1, 경계모수 간격=1.5, 심각도=3, 부적편포

### 문항 1 (조작 문항)

| 항목 | 내용 |
|------|------|
| 채점함수 | `[0, 1, 2, sf4, 4]` — sf4는 IV1 수준에 따라 3 / 3.33 / 3.66 |
| 경계모수 | `[μ−1.5d, μ−0.5d, μ+0.5d, μ+1.5d]` — 평균=μ(IV3), 간격=d(IV2) |
| 변별도 | 1 (고정) |

**경계모수 예시** (μ=0, d=2): `[−3, −1, 1, 3]` → 평균 = 0 ✓

### 문항 2–20 (고정 구조)

| 항목 | 내용 |
|------|------|
| 채점함수 | `[0, 1, 2, 3, 4]` (등간격 고정) |
| 경계모수 간격 | 1.5 고정 |
| 경계모수 평균 | 반복마다 `N(0, 1)` 에서 독립 표집 |
| 변별도 | 1 (고정) |

### 능력모수 분포

| 수준 | 방법 | 평균 | 분산 |
|------|------|------|------|
| 정적편포 | `sn::rsn()`, skewness ≈ +0.8 | ≈ 0 | ≈ 1 |
| 정규분포 | `rnorm(500, 0, 1)` | 0 | 1 |
| 부적편포 | `sn::rsn()`, skewness ≈ −0.8 | ≈ 0 | ≈ 1 |
| 균등분포 | `runif(500, −4, 4)` | 0 | ≈ 5.33 |

> 균등분포는 분산이 다른 분포와 다름 — 의도적 설계

`sn::rsn()`의 `alpha` 파라미터는 왜도 값과 직접 일치하지 않으므로, `00_setup.R`에서 수치 최적화(`uniroot`)로 평균=0, 분산=1을 만족하는 `xi`, `omega`, `alpha`를 사전 계산한다.

### 모수 추정

- 패키지: `mirt`
- 모형: 1차원 GPCM (표준 채점함수 [0,1,2,3,4] 가정)
- 의도: 비등간격 채점함수로 생성된 데이터를 등간격으로 잘못 가정하여 추정

---

## 실행 방법

### 1단계: 시드 테이블 생성 (최초 1회)

```r
setwd("simulation1")   # 작업 디렉토리를 simulation1/로 설정
source("R/01_seed_generator.R")
```

`output/seed_table.rds` 생성 후 이 파일은 **절대 삭제하지 말 것** (재현성의 근거).

### 2단계: 조건별 실행

`R/06_main.R` 상단 **사용자 설정 영역**을 수정 후 실행:

```r
COND_CODE    <- "1233"  # 실행할 조건 코드
USE_PARALLEL <- FALSE   # 병렬 여부
N_CORES      <- 4       # 병렬 코어 수
N_REP        <- 1000    # 반복 횟수
```

```r
source("R/06_main.R")
```

### 다른 컴퓨터에서 이어 실행하기

중단 후 다른 컴퓨터에서 재개할 때 반드시 아래 파일을 옮겨야 한다:

```
output/seed_table.rds                  ← 필수 (시드 재현성)
output/progress/cond[C]_completed.txt  ← 중단된 조건의 진행 상황
```

옮긴 후 동일한 `COND_CODE`로 `06_main.R`을 실행하면, 완료된 반복은 자동으로 건너뛰고 미완료 반복만 이어서 실행된다.

---

## 저장 파일 설명

### `output/true_params/..._item_params.csv`

| 열 | 내용 |
|----|------|
| item | 문항 번호 (1–20) |
| a | 변별도 |
| b1–b4 | 경계 모수 |
| sf1–sf5 | 채점함수 값 |

### `output/responses/..._response.csv`

| 열 | 내용 |
|----|------|
| theta | 능력모수 (500행) |
| item1–item20 | 응답값 (0–4) |

### `output/estimated_params/..._est_params.csv`

성공 시 (20행):

| 열 | 내용 |
|----|------|
| item | 문항명 (Item1 등) |
| converged | 수렴 여부 |
| a | 추정 변별도 |
| b1–b4 | 추정 경계 모수 |

실패 시 (1행): `success=FALSE, error=에러메시지`

### 로그 파일 형식

```
[YYYY-MM-DD HH:MM:SS] rep=NNNN func=FUNC_NAME                   status=OK    msg=...
[2026-06-25 14:03:11] rep=0001 func=generate_item_params         status=OK    msg=
[2026-06-25 14:03:12] rep=0001 func=generate_response            status=OK    msg=
[2026-06-25 14:03:15] rep=0001 func=estimate_params              status=OK    msg=converged=TRUE
```

---

## 유의사항

1. **작업 디렉토리**: `simulation1/` 을 작업 디렉토리로 설정 후 실행할 것  
   (`setwd("경로/simulation1")` 또는 RStudio에서 `simulation1.Rproj` 열기)

2. **시드 테이블 보존**: `output/seed_table.rds` 를 삭제하면 재현성이 깨짐

3. **병렬 실행 시 주의**:
   - 로그는 `output/logs/temp/` 에 임시 저장 후, 실행 완료 시 자동 병합됨
   - 병렬 도중 강제 종료 시 `.tmp` 파일이 남을 수 있음 → 수동으로 병합 필요:
     ```r
     source("R/05_logger.R")
     merge_temp_logs("output/logs/gen_log_cond1111.txt")
     merge_temp_logs("output/logs/est_log_cond1111.txt")
     ```

4. **응답값 범위**: 가상 응답은 0–4 범위로 저장, mirt 추정 시 내부에서 +1 변환

5. **재현성 확인**: 동일 조건·반복을 재실행하면 동일한 시드로 동일 결과 생성됨  
   (단, `output/seed_table.rds` 와 R 버전·패키지 버전이 같아야 함)

6. **패키지 버전 기록**: 분석 보고 시 아래 정보 포함 권장
   ```r
   sessionInfo()
   packageVersion("mirt")
   packageVersion("sn")
   ```
