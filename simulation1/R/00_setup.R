# =============================================================================
# 00_setup.R
# 패키지 로드, 전역 상수 정의, rsn 보정 파라미터 계산
# =============================================================================

# ── 패키지 ──────────────────────────────────────────────────────────────────
required_pkgs <- c("mirt", "sn", "doParallel", "foreach")
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
# 조건코드 첫째 자리: 채점함수 네 번째 값 (sf[4])
IV1_LEVELS <- c(3, 3.33, 3.66)

# 조건코드 둘째 자리: 경계 모수 간격
IV2_LEVELS <- c(2, 1.5, 1)

# 조건코드 셋째 자리: 문항 심각도(경계 모수 평균)
IV3_LEVELS <- c(0, 1.5, 3)

# 조건코드 넷째 자리: 능력모수 분포
IV4_LEVELS <- c("pos_skew", "normal", "neg_skew", "uniform")

# ── rsn 보정 파라미터 계산 (평균=0, 분산=1) ──────────────────────────────────
# sn 패키지 rsn(xi, omega, alpha)에서
#   E[X]   = xi + omega * delta * sqrt(2/pi)
#   Var[X] = omega^2 * (1 - 2*delta^2/pi)
#   delta  = alpha / sqrt(1 + alpha^2)
#
# 목표: E[X]=0, Var[X]=1 이 되도록 xi, omega 역산
# 주의: alpha != 왜도(skewness); alpha는 수치 최적화로 결정

find_rsn_params <- function(target_skew) {
  # Step 1: target_skew를 만족하는 alpha 탐색
  f_alpha <- function(alpha) {
    delta    <- alpha / sqrt(1 + alpha^2)
    b_const  <- sqrt(2 / pi)
    skewness <- ((4 - pi) / 2) * (delta * b_const)^3 /
                (1 - 2 * delta^2 / pi)^(3 / 2)
    skewness - target_skew
  }
  interval <- if (target_skew > 0) c(0.001, 50) else c(-50, -0.001)
  alpha    <- uniroot(f_alpha, interval = interval, tol = 1e-10)$root

  # Step 2: 평균=0, 분산=1이 되도록 xi, omega 역산
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- 1 / sqrt(1 - 2 * delta^2 / pi)
  xi    <- -omega * delta * sqrt(2 / pi)

  list(alpha = alpha, xi = xi, omega = omega)
}

RSN_POS <- find_rsn_params(0.8)   # 정적편포 (skewness ≈ +0.8)
RSN_NEG <- find_rsn_params(-0.8)  # 부적편포 (skewness ≈ -0.8)

# ── 모든 유효 조건코드 목록 ───────────────────────────────────────────────────
ALL_COND_CODES <- character(0)
for (i1 in 1:3) for (i2 in 1:3) for (i3 in 1:3) for (i4 in 1:4) {
  ALL_COND_CODES <- c(ALL_COND_CODES, paste0(i1, i2, i3, i4))
}
# 총 108개

cat("setup.R 로드 완료 — 조건 수:", length(ALL_COND_CODES), "\n")
cat(sprintf("  RSN_POS: alpha=%.4f, xi=%.4f, omega=%.4f\n",
            RSN_POS$alpha, RSN_POS$xi, RSN_POS$omega))
cat(sprintf("  RSN_NEG: alpha=%.4f, xi=%.4f, omega=%.4f\n",
            RSN_NEG$alpha, RSN_NEG$xi, RSN_NEG$omega))
