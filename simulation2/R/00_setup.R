# =============================================================================
# 00_setup.R  (simulation2)
# 패키지 로드, 전역 상수 정의
#
# simulation2 변경 사항:
#   - IV1: 7 수준 (sf4 값 = 3, 3.2, 3.33, 3.4, 3.6, 3.66, 3.8)
#   - IV2: 10 수준 (경계모수 간격 = 1.4, 1.25, 1.2, 1, 0.8, 0.6, 0.5, 0.4, 0.25, 0.2)
#   - IV3: 4 수준 (경계모수 평균 = 0, 0.5, 1, 1.5)
#   - IV4: 1 수준 (정규분포만)
#   - 조건 코드: 8자리 2자리 고정 형식 (예: "01010101")
#   - 총 조건 수: 7 × 10 × 4 × 1 = 280
# =============================================================================

# ── 패키지 ──────────────────────────────────────────────────────────────────
required_pkgs <- c("mirt", "doParallel", "foreach")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
  library(pkg, character.only = TRUE)
}

# ── 출력 디렉토리 생성 ───────────────────────────────────────────────────────
dirs <- c(
  "output/logs/temp",
  "output/progress",
  "output/true_params",
  "output/responses",
  "output/estimated_params"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── 고정 설계 상수 ───────────────────────────────────────────────────────────
N_ITEMS   <- 20
N_CAT     <- 5
N_PERSONS <- 500
DISCRIM   <- 1

# ── 독립변수 수준 매핑 ───────────────────────────────────────────────────────
# IV1: 채점함수 네 번째 값 (sf[4]) — 7 수준
IV1_LEVELS <- c(3, 3.2, 3.33, 3.4, 3.6, 3.66, 3.8)

# IV2: 경계 모수 간격 — 10 수준
IV2_LEVELS <- c(1.4, 1.25, 1.2, 1, 0.8, 0.6, 0.5, 0.4, 0.25, 0.2)

# IV3: 문항 심각도 (경계 모수 평균) — 4 수준
IV3_LEVELS <- c(0, 0.5, 1, 1.5)

# IV4: 능력모수 분포 — 1 수준 (정규분포만)
IV4_LEVELS <- c("normal")

# ── 조건 코드 생성 (2자리 고정 형식) ─────────────────────────────────────────
# 형식: "AABBCCDD" — AA=IV1 인덱스, BB=IV2 인덱스, CC=IV3 인덱스, DD=IV4 인덱스
# 예: IV1=1번째, IV2=3번째, IV3=2번째, IV4=1번째 → "01030201"
ALL_COND_CODES <- character(0)
for (i1 in seq_along(IV1_LEVELS))
  for (i2 in seq_along(IV2_LEVELS))
    for (i3 in seq_along(IV3_LEVELS))
      for (i4 in seq_along(IV4_LEVELS))
        ALL_COND_CODES <- c(ALL_COND_CODES,
                            sprintf("%02d%02d%02d%02d", i1, i2, i3, i4))
# 총 7 × 10 × 4 × 1 = 280개

cat("setup.R 로드 완료 — 조건 수:", length(ALL_COND_CODES), "\n")
cat(sprintf("  IV1 (%d 수준): %s\n", length(IV1_LEVELS),
            paste(IV1_LEVELS, collapse = ", ")))
cat(sprintf("  IV2 (%d 수준): %s\n", length(IV2_LEVELS),
            paste(IV2_LEVELS, collapse = ", ")))
cat(sprintf("  IV3 (%d 수준): %s\n", length(IV3_LEVELS),
            paste(IV3_LEVELS, collapse = ", ")))
cat(sprintf("  IV4 (%d 수준): %s\n", length(IV4_LEVELS),
            paste(IV4_LEVELS, collapse = ", ")))
