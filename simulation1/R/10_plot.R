# =============================================================================
# 10_plot.R
# 기존 분석 결과 CSV를 읽어 그래프를 재생성하는 스크립트
#
# 입력:
#   output/analysis/transposition_summary.csv
#   output/analysis/rmse_bias_summary.csv
#
# 출력:
#   output/analysis/transposition_plot.png
#   output/analysis/rmse_bias_plot_b_bias.png
#   output/analysis/rmse_bias_plot_b_rmse.png
#
# 변경 사항 (기존 분석 파일 대비):
#   - 모든 텍스트 영어
#   - 제목 제거
#   - 글씨 크기 확대 (base_size = 14)
#   - 채점함수 약어: sf → s
#   - b-parameter Bias/RMSE: facet 열 순서 IV1↔IV2 교환
#     (열: Interval 0.5 / 1.0 / 1.5, 그 안에 s4 3.00 / 3.33 / 3.66)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── IV4 영어 라벨 매핑 ────────────────────────────────────────────────────────
IV4_EN_LEVELS <- c("normal", "pos_skew", "neg_skew", "uniform")
IV4_EN_LABELS <- c("Normal", "Negatively Skewed", "Positively Skewed", "Uniform")

# ── factor 헬퍼 ───────────────────────────────────────────────────────────────
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

# =============================================================================
# 1. transposition_plot
# =============================================================================
trans_path <- "output/analysis/transposition_summary.csv"
if (!file.exists(trans_path)) {
  message("건너뜀: ", trans_path, " 없음")
} else {
  df_t <- read.csv(trans_path, stringsAsFactors = FALSE)

  df_t <- df_t %>%
    mutate(
      iv1 = make_num_factor(iv1_sf4,        "%.2f", decreasing = FALSE),
      iv2 = make_num_factor(iv2_b_interval, "%.1f", decreasing = FALSE),
      iv3 = make_num_factor(iv3_b_mean,     "%.1f", decreasing = FALSE),
      iv4 = factor(iv4_theta_dist,
                   levels = IV4_EN_LEVELS,
                   labels = IV4_EN_LABELS)
    )

  p_trans <- ggplot(df_t,
                    aes(x    = iv3,
                        y    = prop_transposed,
                        fill = iv4)) +
    geom_col(position = "dodge") +
    facet_grid(
      rows     = vars(iv1),
      cols     = vars(iv2),
      labeller = labeller(
        iv1 = function(x) paste0("s4 = ", x),
        iv2 = function(x) paste0("Interval = ", x)
      )
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x    = "b-parameter Mean (IV3)",
      y    = "Reversal Rate",
      fill = "Ability Distribution (IV4)"
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(size = 13),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 14),
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 13)
    )

  ggsave("output/analysis/transposition_plot.png", p_trans,
         width = 10, height = 8, dpi = 150)
  message("저장 완료: output/analysis/transposition_plot.png")
}

# =============================================================================
# 2. rmse_bias plots (b1–b4)
# =============================================================================
rmse_path <- "output/analysis/rmse_bias_summary.csv"
if (!file.exists(rmse_path)) {
  message("건너뜀: ", rmse_path, " 없음")
} else {
  df_r <- read.csv(rmse_path, stringsAsFactors = FALSE)

  b_params <- c("b1", "b2", "b3", "b4")
  b_labels <- c(b1 = "b₁", b2 = "b₂", b3 = "b₃", b4 = "b₄")

  # ── IV1: iv2를 외곽 열(0.5 → 1.0 → 1.5 왼쪽부터), iv1을 내측 열 ─────────
  # iv2: decreasing=FALSE → 0.5, 1.0, 1.5 순
  # iv1: decreasing=FALSE → 3.00, 3.33, 3.66 순
  plot_b <- df_r %>%
    filter(param %in% b_params) %>%
    mutate(
      param = factor(param, levels = b_params, labels = b_labels),
      iv1   = make_num_factor(iv1_sf4,        "%.2f", decreasing = FALSE),
      iv2   = make_num_factor(iv2_b_interval, "%.1f", decreasing = FALSE),
      iv3   = make_num_factor(iv3_b_mean,     "%.1f", decreasing = FALSE),
      iv4   = factor(iv4_theta_dist,
                     levels = IV4_EN_LEVELS,
                     labels = IV4_EN_LABELS)
    )

  common_layers <- list(
    facet_grid(
      rows     = vars(param),
      cols     = vars(iv2, iv1),
      labeller = labeller(
        iv2   = function(x) paste0("Interval = ", x),
        iv1   = function(x) paste0("s4 = ", x),
        param = label_value
      )
    ),
    scale_colour_brewer(palette = "Set1"),
    labs(
      x      = "b-parameter Mean (IV3)",
      colour = "Ability Distribution (IV4)"
    ),
    theme_bw(base_size = 14),
    theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 11),
      axis.text       = element_text(size = 11),
      axis.title      = element_text(size = 13),
      legend.text     = element_text(size = 11),
      legend.title    = element_text(size = 12)
    )
  )

  # ── Bias 그래프 ──────────────────────────────────────────────────────────────
  p_bias <- ggplot(plot_b,
                   aes(x = iv3, y = bias,
                       colour = iv4, group = iv4)) +
    geom_hline(yintercept = 0,
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    labs(y = "Bias") +
    common_layers

  ggsave("output/analysis/rmse_bias_plot_b_bias.png", p_bias,
         width = 16, height = 10, dpi = 150)
  message("저장 완료: output/analysis/rmse_bias_plot_b_bias.png")

  # ── RMSE 그래프 ──────────────────────────────────────────────────────────────
  p_rmse <- ggplot(plot_b,
                   aes(x = iv3, y = rmse,
                       colour = iv4, group = iv4)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    labs(y = "RMSE") +
    common_layers

  ggsave("output/analysis/rmse_bias_plot_b_rmse.png", p_rmse,
         width = 16, height = 10, dpi = 150)
  message("저장 완료: output/analysis/rmse_bias_plot_b_rmse.png")
}

message("\n10_plot.R 완료")
