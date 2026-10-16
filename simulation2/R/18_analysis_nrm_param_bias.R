# =============================================================================
# 18_analysis_nrm_param_bias.R  (simulation2, NRM 추정 결과)
# NRM 문항 1의 기울기 α, 채점함수 s_k, 절편 c_k 추정의 Bias / RMSE 계산·시각화
#
# ── 모수와 참값 ──────────────────────────────────────────────────────────────
#   mirt NRM :  P(X=k|θ) ∝ exp(α · s_k · θ + c_k),  k = 0..4
#               α = a1,  s_k = ak_k,  c_k = d_k
#               식별 제약(mirt 기본): s_0 = 0, s_4 = 4, c_0 = 0  → 추정 모수: α, s1~s3, c1~c4
#   생성 모형:  η_k = Σ_{j≤k} a · w_j · (θ − b_j),  w_j = sf_j − sf_{j−1},  sf = (0, 1, 2, s4, 4), a = 1
#               = a · sf_k · θ − a · Σ_{j≤k} w_j · b_j
#   → 참값:   α = a = 1
#             s_k = sf_k             (s1 = 1, s2 = 2, s3 = IV1 의 s4 값)
#             c_k = −Σ_{j≤k} w_j b_j  (b_j = 평균 + 간격 × (−1.5, −0.5, 0.5, 1.5))
#   NRM 은 생성 모형을 그대로 포함하므로(오지정 없음) 기준값은 참값 자체
#
# ── 지표 (조건 × 모수별) ──────────────────────────────────────────────────────
#   bias = mean(추정 − 참값),  rmse = sqrt(mean((추정 − 참값)²))
#   median_bias = median(추정 − 참값),  mdae = median(|추정 − 참값|)   ← 극단값에 덜 민감
#
# ── 반복 집합 (16_plot_nrm_excluded.R 과 같은 정의, 반복 단위 제외) ──────────
#   all                — 수렴하고 b1~b4 가 계산된 반복 전체
#   excl_extreme       — b1~b4 중 하나라도 |b̂ − b| ≥ ERR_CUT 인 반복 제외
#   excl_extreme_sfrev — 위 + Δs_k ≤ 0 또는 |Δs_k| < DAK_CUT 인 반복 제외
#     ※ excl_extreme_sfrev 는 s_k 추정값으로 반복을 고르므로 s_k 의 Bias 가 그만큼 작아지는 쪽으로 치우침
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력)
# 출력: output/nrm/param_bias/<반복 집합>/
#   param_bias_summary.csv   — 조건 × 모수별 n, 참값, 추정 평균, Bias, RMSE, 중앙값 편향, MdAE
#                              (Bias·RMSE 는 모수 하나씩 따로 계산 — 모수끼리 합치지 않음)
#   <모수>_bias.png / <모수>_rmse.png — 모수마다 따로 그린 그래프 (모수 8개 × 2 = 16개)
#     alpha_bias.png, alpha_rmse.png              — 기울기 α
#     s1_bias.png ~ s3_bias.png, s1_rmse.png ~ s3_rmse.png — 채점함수 s1, s2, s3
#     c1_bias.png ~ c4_bias.png, c1_rmse.png ~ c4_rmse.png — 절편 c1 ~ c4
#   (배치: 행 IV1(s4), 열 IV2(간격, 좁음→넓음), 패널 x축 IV3(평균) — 전치 그래프와 같음)
#   대상: 문항 1 (조작 문항) 의 모수만 — 11번이 문항 1 모수만 저장하므로 nrm_test_rep.csv 에는 문항 1만 있음
# =============================================================================

ERR_CUT <- 3     # 극단 오차: |b̂ − b| ≥ ERR_CUT   (15·16번과 같은 값)
DAK_CUT <- 0.1   # 분모 붕괴: |Δs_k| < DAK_CUT      (15·16번과 같은 값)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
BASE_DIR <- "output/nrm/param_bias"

B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

# =============================================================================
# 1. 데이터와 참값
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

# 참 경계모수와 참 NRM 모수 (반복마다 같은 값이지만 조건별로 계산)
bt <- sapply(B_OFFSET, function(o) rep_df$iv3_b_mean + rep_df$iv2_b_interval * o)  # n × 4
sf <- cbind(0, 1, 2, rep_df$iv1_sf4, 4)
w  <- sf[, 2:5] - sf[, 1:4]                                                         # n × 4
c_true <- t(apply(w * bt, 1, cumsum)) * -1                                          # c1~c4

rep_df <- rep_df %>%
  mutate(
    alpha_true = 1,
    s1_true = 1, s2_true = 2, s3_true = iv1_sf4,
    c1_true = c_true[, 1], c2_true = c_true[, 2], c3_true = c_true[, 3], c4_true = c_true[, 4],
    alpha_est = a1,
    s1_est = ak1, s2_est = ak2, s3_est = ak3,
    c1_est = d1, c2_est = d2, c3_est = d3, c4_est = d4,
    # 16번과 같은 제외 기준
    is_extreme = pmax(abs(b1 - bt[, 1]), abs(b2 - bt[, 2]),
                      abs(b3 - bt[, 3]), abs(b4 - bt[, 4])) >= ERR_CUT,
    is_sfrev   = {
      dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
      rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
    }
  )

PARAMS <- c("alpha", "s1", "s2", "s3", "c1", "c2", "c3", "c4")

long_all <- rep_df %>%
  select(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean, is_extreme, is_sfrev,
         ends_with("_true"), ends_with("_est")) %>%
  select(-any_of(c("b3_true", "b4_true"))) %>%
  pivot_longer(c(ends_with("_true"), ends_with("_est")),
               names_to = c("param", ".value"), names_sep = "_") %>%
  filter(param %in% PARAMS) %>%
  mutate(err = est - true)

# =============================================================================
# 2. 그래프 도우미
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
                iv2 = make_num_factor(iv2_b_interval, "%.2f"),   # 좁음 → 넓음
                iv3 = make_num_factor(iv3_b_mean,     "%.1f"))
}
PARAM_TITLES <- c(alpha = "alpha (slope, true = 1)",
                  s1 = "s1 (scoring function, true = 1)",
                  s2 = "s2 (scoring function, true = 2)",
                  s3 = "s3 (scoring function, true = s4)",
                  c1 = "c1 (intercept)", c2 = "c2 (intercept)",
                  c3 = "c3 (intercept)", c4 = "c4 (intercept)")

# 모수 하나의 Bias 또는 RMSE 그래프
param_plot <- function(df, pp, yvar, desc, path) {
  pd <- df %>% filter(param == pp)
  is_bias <- yvar == "bias"
  colour  <- if (is_bias) "steelblue" else "firebrick"
  p <- ggplot(pd, aes(x = iv3, y = .data[[yvar]], group = 1))
  if (is_bias) p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4)
  p <- p + geom_line(linewidth = 0.7, colour = colour) + geom_point(size = 1.4, colour = colour) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    labs(title = sprintf("Item 1 NRM: %s - %s", PARAM_TITLES[[pp]], desc),
         x = "b-parameter Mean (IV3)",
         y = sprintf("%s of %s", if (is_bias) "Bias" else "RMSE", pp)) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path, p, width = 24, height = 18, dpi = 150)
}

# =============================================================================
# 3. 반복 집합별 계산·그래프
# =============================================================================
run_set <- function(keep_fun, label, desc) {
  out_dir <- file.path(BASE_DIR, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)

  ld <- long_all %>% filter(keep_fun(is_extreme, is_sfrev), is.finite(err))
  summary_df <- ld %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
    summarise(n           = n(),
              true        = first(true),
              mean_est    = mean(est),
              bias        = mean(err),
              rmse        = sqrt(mean(err^2)),
              median_bias = median(err),
              mdae        = median(abs(err)),
              .groups = "drop") %>%
    mutate(param = factor(param, levels = PARAMS)) %>%
    arrange(cond_code, param)
  write.csv(summary_df, out("param_bias_summary.csv"), row.names = FALSE)

  pdf <- add_factors(summary_df %>% mutate(param = as.character(param)))

  # 모수마다 따로: Bias 1장, RMSE 1장
  for (pp in PARAMS) {
    param_plot(pdf, pp, "bias", desc, out(sprintf("%s_bias.png", pp)))
    param_plot(pdf, pp, "rmse", desc, out(sprintf("%s_rmse.png", pp)))
  }

  n_rep <- n_distinct(paste(ld$cond_code, ld$rep_id))
  message(sprintf("[%s] 반복 %d개 → %s/", label, n_rep, out_dir))
}

run_set(function(ext, sfr) rep(TRUE, length(ext)), "all",
        "all converged replications")
run_set(function(ext, sfr) !ext, "excl_extreme",
        sprintf("excluding extreme errors (|error| >= %g)", ERR_CUT))
run_set(function(ext, sfr) !ext & !sfr, "excl_extreme_sfrev",
        sprintf("excluding extreme errors and SF reversals / near-zero ds (|ds| < %g)", DAK_CUT))
