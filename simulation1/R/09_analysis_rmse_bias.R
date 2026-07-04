# =============================================================================
# 08_analysis_rmse_bias.R
# 문항 1 추정 모수의 편향(Bias)과 RMSE 분석
#
# 분석 대상: 문항 1 (조작 문항) — a, b1, b2, b3, b4
#
# 참값 계산:
#   a_true  = 1 (DISCRIM 고정)
#   b_true  = iv3_b_mean + iv2_b_interval × c(-1.5, -0.5, 0.5, 1.5)
#   (반복마다 동일한 결정론적 값 — 랜덤 성분 없음)
#
# 출력 파일:
#   output/analysis/rmse_bias_long.csv       — 반복별 오차 원자료 (long format)
#   output/analysis/rmse_bias_summary.csv    — 조건 × 모수별 Bias / RMSE 집계
#   output/analysis/rmse_bias_plot_b.png     — b1-b4 Bias / RMSE 시각화
#   output/analysis/rmse_bias_plot_a.png     — a 모수 Bias / RMSE 시각화
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── 조건 코드 → 실제 값 매핑 (시뮬레이션 설계 변경 시 여기만 수정) ──────────
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")
IV4_LABEL <- c(
  "pos_skew" = "정적편포(+0.8)",
  "normal"   = "정규분포",
  "neg_skew" = "부적편포(-0.8)",
  "uniform"  = "균등분포"
)
DISCRIM <- 1   # 고정 변별도

# =============================================================================
# 1단계: 추정 모수 CSV 로딩 — 문항 1 행만 추출
# =============================================================================
est_files <- list.files(
  path       = "output/estimated_params",
  pattern    = "_est_params\\.csv$",
  full.names = TRUE
)
if (length(est_files) == 0) stop("추정 결과 파일이 없습니다. 먼저 시뮬레이션을 실행하세요.")
cat(sprintf("파일 %d개 로딩 중...\n", length(est_files)))

parse_fname <- function(path) {
  fname <- basename(path)
  cond  <- str_match(fname, "_cond([0-9]+)_")[, 2]
  rep   <- as.integer(str_match(fname, "_rep([0-9]+)_")[, 2])
  list(cond_code = cond, rep_id = rep)
}

read_item1_est <- function(path) {
  meta <- parse_fname(path)
  df   <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(df)) return(NULL)

  # 추정 실패 파일
  if ("success" %in% names(df) && isFALSE(df$success[1])) {
    return(data.frame(
      cond_code = meta$cond_code, rep_id = meta$rep_id,
      converged = FALSE,
      a_est = NA_real_,
      b1_est = NA_real_, b2_est = NA_real_,
      b3_est = NA_real_, b4_est = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # 문항 1 행 추출
  item1 <- df[df$item %in% c("Item1", "item1", 1), ]
  if (nrow(item1) == 0) item1 <- df[1, ]

  data.frame(
    cond_code = meta$cond_code,
    rep_id    = meta$rep_id,
    converged = isTRUE(item1$converged[1]),
    a_est     = item1$a[1],
    b1_est    = item1$b1[1],
    b2_est    = item1$b2[1],
    b3_est    = item1$b3[1],
    b4_est    = item1$b4[1],
    stringsAsFactors = FALSE
  )
}

est_list <- lapply(est_files, read_item1_est)
est_raw  <- do.call(rbind, est_list[!sapply(est_list, is.null)])
cat(sprintf("  총 %d행 로딩 완료\n", nrow(est_raw)))

# =============================================================================
# 2단계: 조건 코드 분해 및 참값 계산
# =============================================================================
est_raw <- est_raw %>%
  mutate(
    d1 = substr(cond_code, 1, 1),
    d2 = substr(cond_code, 2, 2),
    d3 = substr(cond_code, 3, 3),
    d4 = substr(cond_code, 4, 4),
    iv1_sf4        = IV1_MAP[d1],
    iv2_b_interval = IV2_MAP[d2],
    iv3_b_mean     = IV3_MAP[d3],
    iv4_theta_dist = IV4_MAP[d4],
    iv4_label      = IV4_LABEL[IV4_MAP[d4]],
    # 문항 1 참값 (결정론적)
    a_true  = DISCRIM,
    b1_true = iv3_b_mean + iv2_b_interval * (-1.5),
    b2_true = iv3_b_mean + iv2_b_interval * (-0.5),
    b3_true = iv3_b_mean + iv2_b_interval *   0.5,
    b4_true = iv3_b_mean + iv2_b_interval *   1.5
  ) %>%
  select(-d1, -d2, -d3, -d4)

# =============================================================================
# 3단계: 오차 계산 (추정 - 참값)
# =============================================================================
est_raw <- est_raw %>%
  mutate(
    err_a  = a_est  - a_true,
    err_b1 = b1_est - b1_true,
    err_b2 = b2_est - b2_true,
    err_b3 = b3_est - b3_true,
    err_b4 = b4_est - b4_true
  )

# =============================================================================
# 4단계: Long format 오차 원자료 저장
# =============================================================================
long_df <- est_raw %>%
  filter(converged == TRUE) %>%
  select(cond_code, rep_id,
         iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist, iv4_label,
         a_est, b1_est, b2_est, b3_est, b4_est,
         a_true, b1_true, b2_true, b3_true, b4_true,
         err_a, err_b1, err_b2, err_b3, err_b4) %>%
  arrange(cond_code, rep_id)

write.csv(long_df, "output/analysis/rmse_bias_long.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/rmse_bias_long.csv (%d행)\n", nrow(long_df)))

# =============================================================================
# 5단계: 조건 × 모수별 Bias / RMSE 집계
# =============================================================================
group_vars <- c("cond_code", "iv1_sf4", "iv2_b_interval",
                "iv3_b_mean", "iv4_theta_dist", "iv4_label")

compute_metrics <- function(df, err_col, param_name) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n          = n(),
      bias       = mean(.data[[err_col]], na.rm = TRUE),
      rmse       = sqrt(mean(.data[[err_col]]^2, na.rm = TRUE)),
      sd_err     = sd(.data[[err_col]], na.rm = TRUE),
      .groups    = "drop"
    ) %>%
    mutate(param = param_name)
}

summary_df <- bind_rows(
  compute_metrics(long_df, "err_a",  "a"),
  compute_metrics(long_df, "err_b1", "b1"),
  compute_metrics(long_df, "err_b2", "b2"),
  compute_metrics(long_df, "err_b3", "b3"),
  compute_metrics(long_df, "err_b4", "b4")
) %>%
  select(cond_code, param, everything()) %>%
  arrange(cond_code, param)

write.csv(summary_df, "output/analysis/rmse_bias_summary.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/rmse_bias_summary.csv (%d행)\n", nrow(summary_df)))

# =============================================================================
# 6단계: 시각화 공통 설정
# =============================================================================
make_iv_factor <- function(x, fmt = "%.2f", decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

summary_df <- summary_df %>%
  mutate(
    iv1 = make_iv_factor(iv1_sf4,        fmt = "%.2f", decreasing = FALSE),
    iv2 = make_iv_factor(iv2_b_interval, fmt = "%.1f", decreasing = TRUE),
    iv3 = make_iv_factor(iv3_b_mean,     fmt = "%.1f", decreasing = FALSE),
    iv4 = factor(iv4_theta_dist,
                 levels = c("normal", "pos_skew", "neg_skew", "uniform"),
                 labels = c("정규분포", "정적편포(+)", "부적편포(-)", "균등분포")),
    iv4_label = factor(iv4_label,
                       levels = c("정규분포", "정적편포(+0.8)",
                                  "부적편포(-0.8)", "균등분포"))
  )

# =============================================================================
# 7단계: b1–b4 Bias / RMSE 시각화
# =============================================================================
b_params <- c("b1", "b2", "b3", "b4")
b_labels <- c(b1 = "b₁", b2 = "b₂", b3 = "b₃", b4 = "b₄")

plot_b <- summary_df %>%
  filter(param %in% b_params) %>%
  mutate(param = factor(param, levels = b_params, labels = b_labels))

# ── Bias 그래프 ──
p_bias_b <- ggplot(plot_b,
                   aes(x = iv3, y = bias, colour = iv4_label, group = iv4_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(
    rows = vars(param),
    cols = vars(iv1, iv2),
    labeller = labeller(
      iv1   = function(x) paste0("sf4=", x),
      iv2   = function(x) paste0("간격=", x),
      param = label_value
    )
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title   = "문항 1 경계 모수 Bias (추정 - 참값)",
    x       = "문항 심각도 — IV3: 경계모수 평균",
    y       = "Bias",
    colour  = "능력모수 분포 (IV4)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text      = element_text(size = 8))

ggsave("output/analysis/rmse_bias_plot_b_bias.png", p_bias_b,
       width = 14, height = 10, dpi = 150)
cat("저장 완료: output/analysis/rmse_bias_plot_b_bias.png\n")

# ── RMSE 그래프 ──
p_rmse_b <- ggplot(plot_b,
                   aes(x = iv3, y = rmse, colour = iv4_label, group = iv4_label)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(
    rows = vars(param),
    cols = vars(iv1, iv2),
    labeller = labeller(
      iv1   = function(x) paste0("sf4=", x),
      iv2   = function(x) paste0("간격=", x),
      param = label_value
    )
  ) +
  scale_colour_brewer(palette = "Set1") +
  labs(
    title   = "문항 1 경계 모수 RMSE",
    x       = "문항 심각도 — IV3: 경계모수 평균",
    y       = "RMSE",
    colour  = "능력모수 분포 (IV4)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text      = element_text(size = 8))

ggsave("output/analysis/rmse_bias_plot_b_rmse.png", p_rmse_b,
       width = 14, height = 10, dpi = 150)
cat("저장 완료: output/analysis/rmse_bias_plot_b_rmse.png\n")

# =============================================================================
# 8단계: a 모수 Bias / RMSE 시각화
# =============================================================================
plot_a <- summary_df %>% filter(param == "a")

p_a <- ggplot(plot_a,
              aes(x = iv3, colour = iv4_label, group = iv4_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(aes(y = bias),  linewidth = 0.7) +
  geom_point(aes(y = bias), size = 1.8, shape = 16) +
  geom_line(aes(y = rmse),  linewidth = 0.7, linetype = "dotted") +
  geom_point(aes(y = rmse), size = 1.8, shape = 17) +
  facet_grid(
    rows = vars(iv2),
    cols = vars(iv1),
    labeller = labeller(
      iv1 = function(x) paste0("sf4=", x),
      iv2 = function(x) paste0("간격=", x)
    )
  ) +
  scale_colour_brewer(palette = "Set1") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = "실선=Bias  점선=RMSE", size = 2.8, colour = "grey40") +
  labs(
    title   = "문항 1 변별도(a) Bias / RMSE",
    x       = "문항 심각도 — IV3: 경계모수 평균",
    y       = "값",
    colour  = "능력모수 분포 (IV4)"
  ) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text      = element_text(size = 8))

ggsave("output/analysis/rmse_bias_plot_a.png", p_a,
       width = 10, height = 8, dpi = 150)
cat("저장 완료: output/analysis/rmse_bias_plot_a.png\n")

# =============================================================================
# 최종 요약
# =============================================================================
cat("\n================================================================\n")
cat("전체 분석 완료\n")
cat("================================================================\n")
cat(sprintf("  분석 반복 수 (수렴 성공): %d\n", nrow(long_df)))
cat(sprintf("  조건 수                 : %d\n", n_distinct(summary_df$cond_code)))

for (p in c("a", "b1", "b2", "b3", "b4")) {
  sub <- summary_df[summary_df$param == p, ]
  cat(sprintf("  [%s]  Bias 범위: %+.4f ~ %+.4f  |  RMSE 범위: %.4f ~ %.4f\n",
              p,
              min(sub$bias, na.rm = TRUE), max(sub$bias, na.rm = TRUE),
              min(sub$rmse, na.rm = TRUE), max(sub$rmse, na.rm = TRUE)))
}
cat("\n출력 파일:\n")
cat("  output/analysis/rmse_bias_long.csv\n")
cat("  output/analysis/rmse_bias_summary.csv\n")
cat("  output/analysis/rmse_bias_plot_b_bias.png\n")
cat("  output/analysis/rmse_bias_plot_b_rmse.png\n")
cat("  output/analysis/rmse_bias_plot_a.png\n")
