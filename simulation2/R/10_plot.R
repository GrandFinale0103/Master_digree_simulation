# =============================================================================
# 10_plot.R  (simulation2)
# 기존 분석 결과 CSV를 읽어 그래프를 재생성
#
# simulation2 변경:
#   - IV1: 7 수준 (3, 3.2, 3.33, 3.4, 3.6, 3.66, 3.8)
#   - IV2: 10 수준 (1.4 ~ 0.2)
#   - IV3: 4 수준 (0, 0.5, 1, 1.5)
#   - IV4: 정규분포만 (fill/colour 불필요 → 단일 색상)
#
# 입력:
#   output/analysis/transposition_summary.csv
#   output/analysis/rmse_bias_summary.csv
#
# 출력:
#   output/analysis/transposition_plot.png
#   output/analysis/rmse_bias_plot_b_bias.png
#   output/analysis/rmse_bias_plot_b_rmse.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# factor 헬퍼
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
      iv1 = make_num_factor(iv1_sf4,        "%.4f", decreasing = FALSE),
      iv2 = make_num_factor(iv2_b_interval, "%.4f", decreasing = FALSE),
      iv3 = make_num_factor(iv3_b_mean,     "%.2f", decreasing = FALSE)
    )

  # IV4는 정규분포만 → fill 불필요 (단일 색상 막대)
  p_trans <- ggplot(df_t,
                    aes(x = iv3, y = prop_transposed)) +
    geom_col(fill = "steelblue") +
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
    labs(
      x = "b-parameter Mean (IV3)",
      y = "Reversal Rate"
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.text  = element_text(size = 8),
      axis.text   = element_text(size = 8),
      axis.title  = element_text(size = 10)
    )

  ggsave("output/analysis/transposition_plot.png", p_trans,
         width = 20, height = 14, dpi = 150)
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

  # IV2 외곽(좁은 간격→오른쪽), IV1 내측
  # decreasing=TRUE: 간격 큰 것(1.4)이 왼쪽, 작은 것(0.2)이 오른쪽
  plot_b <- df_r %>%
    filter(param %in% b_params) %>%
    mutate(
      param = factor(param, levels = b_params, labels = b_labels),
      iv1   = make_num_factor(iv1_sf4,        "%.4f", decreasing = FALSE),
      iv2   = make_num_factor(iv2_b_interval, "%.4f", decreasing = TRUE),
      iv3   = make_num_factor(iv3_b_mean,     "%.2f", decreasing = FALSE)
    )

  common_layers <- list(
    facet_grid(
      rows     = vars(param),
      cols     = vars(iv2, iv1),
      labeller = labeller(
        iv2   = function(x) paste0("Int=", x),
        iv1   = function(x) paste0("s4=", x),
        param = label_value
      )
    ),
    labs(x = "b-parameter Mean (IV3)"),
    theme_bw(base_size = 9),
    theme(
      strip.text  = element_text(size = 6),
      axis.text   = element_text(size = 7),
      axis.title  = element_text(size = 9)
    )
  )

  # Bias
  p_bias <- ggplot(plot_b, aes(x = iv3, y = bias, group = 1)) +
    geom_hline(yintercept = 0,
               linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7, colour = "steelblue") +
    geom_point(size = 1.5, colour = "steelblue") +
    labs(y = "Bias") +
    common_layers

  ggsave("output/analysis/rmse_bias_plot_b_bias.png", p_bias,
         width = 28, height = 10, dpi = 150)
  message("저장 완료: output/analysis/rmse_bias_plot_b_bias.png")

  # RMSE
  p_rmse <- ggplot(plot_b, aes(x = iv3, y = rmse, group = 1)) +
    geom_line(linewidth = 0.7, colour = "firebrick") +
    geom_point(size = 1.5, colour = "firebrick") +
    labs(y = "RMSE") +
    common_layers

  ggsave("output/analysis/rmse_bias_plot_b_rmse.png", p_rmse,
         width = 28, height = 10, dpi = 150)
  message("저장 완료: output/analysis/rmse_bias_plot_b_rmse.png")
}

message("\n10_plot.R 완료")
