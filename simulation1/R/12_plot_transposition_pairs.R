# =============================================================================
# 12_plot_transposition_pairs.R  (simulation1, GPCM 추정 결과)
# 전치가 어느 경계모수 쌍(b1&b2, b2&b3, b3&b4)에서 발생했는지 시각화
#
# 배치: 원래 전치 비율 그래프(08_analysis_transposition.R)와 동일
#   행 = IV1(s4), 열 = IV2(간격, 왼쪽 좁음 → 오른쪽 넓음), 패널 x축 = IV3(평균)
#   원래 그래프에서 막대 색으로 쓰던 IV4 자리에 경계모수 쌍이 들어가므로
#   IV4(능력모수 분포)는 분포마다 그림을 따로 저장
#
# 입력: output/analysis/transposition_summary.csv  (08_analysis_transposition.R 출력)
#   n_12 / n_23 / n_34 : 해당 쌍이 전치된 반복 수 (한 반복에서 여러 쌍 중복 가능)
#   비율 분모          : 수렴한 반복 수 (n_converged)
#
# 출력:
#   output/analysis/transposition_pairs_summary.csv           — 조건별 쌍별 전치 수·비율
#   output/analysis/transposition_pairs_bar_<분포>.png        — 분포별 막대 그래프
#   output/analysis/transposition_pairs_line_<분포>.png       — 분포별 선 그래프
#   output/analysis/transposition_pairs_stack_<분포>.png      — 분포별 누적 막대 그래프
#     <분포> = normal / pos_skew / neg_skew / uniform
#   output/analysis_normal/transposition_pairs_normal_bar.png — 정규분포 막대 그래프
#   output/analysis_normal/transposition_pairs_normal.png     — 정규분포 선 그래프
#   output/analysis_normal/transposition_pairs_normal_stack.png — 정규분포 누적 막대 그래프
#
# 누적 막대: 쌍별 비율을 쌓은 값이므로, 한 반복에서 여러 쌍이 동시에 전치되면
#            합이 전체 전치 비율보다 커짐 → 전체 전치 비율(prop_transposed)을 검은 점으로 표시
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/analysis/transposition_summary.csv"
if (!file.exists(in_path)) stop("먼저 R/08_analysis_transposition.R 을 실행하세요: ", in_path)

make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

IV4_LEVELS <- c("normal", "pos_skew", "neg_skew", "uniform")
IV4_LABELS <- c("Normal", "Positively Skewed (+0.8)", "Negatively Skewed (-0.8)", "Uniform")

PAIR_LEVELS  <- c("prop_12", "prop_23", "prop_34")
PAIR_LABELS  <- c("b1 & b2", "b2 & b3", "b3 & b4")
pair_colours <- c("b1 & b2" = "#4C9BE8", "b2 & b3" = "#F2A541", "b3 & b4" = "#D1495B")

pairs_df <- read.csv(in_path, stringsAsFactors = FALSE) %>%
  mutate(
    prop_12 = n_12 / n_converged,
    prop_23 = n_23 / n_converged,
    prop_34 = n_34 / n_converged
  )

write.csv(pairs_df %>%
            select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
                   n_converged, n_transposed, prop_transposed,
                   n_12, n_23, n_34, prop_12, prop_23, prop_34),
          "output/analysis/transposition_pairs_summary.csv", row.names = FALSE)

long_df <- pairs_df %>%
  pivot_longer(all_of(PAIR_LEVELS), names_to = "pair", values_to = "prop") %>%
  mutate(
    pair = factor(pair, levels = PAIR_LEVELS, labels = PAIR_LABELS),
    iv1  = make_num_factor(iv1_sf4,        "%.2f"),
    iv2  = make_num_factor(iv2_b_interval, "%.1f", decreasing = FALSE),  # 좁은 → 넓은
    iv3  = make_num_factor(iv3_b_mean,     "%.1f")
  )

# ── 공통 배치: 행 IV1 × 열 IV2, x축 IV3 ──────────────────────────────────────
base_layers <- function(title) list(
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4 = ", x),
                                 iv2 = function(x) paste0("Interval = ", x))),
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)),
  labs(title = title, x = "b-parameter Mean (IV3)", y = "Reversal Rate"),
  theme_bw(base_size = 14),
  theme(legend.position = "bottom", strip.text = element_text(size = 12))
)

plot_bar <- function(df, title) {
  ggplot(df, aes(x = iv3, y = prop, fill = pair)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    scale_fill_manual(values = pair_colours, name = "Reversed pair") +
    base_layers(title)
}
plot_line <- function(df, title) {
  ggplot(df, aes(x = iv3, y = prop, colour = pair, shape = pair, group = pair)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    scale_colour_manual(values = pair_colours, name = "Reversed pair") +
    scale_shape_discrete(name = "Reversed pair") +
    base_layers(title)
}

plot_stack <- function(df, title) {
  total_df <- df %>% distinct(iv1, iv2, iv3, prop_transposed)
  ggplot(df, aes(x = iv3, y = prop)) +
    geom_col(aes(fill = pair), position = position_stack(reverse = TRUE), width = 0.7) +
    geom_point(data = total_df, aes(y = prop_transposed, shape = "Any reversal"),
               size = 2.5, colour = "black") +
    scale_fill_manual(values = pair_colours, name = "Reversed pair") +
    scale_shape_manual(values = c("Any reversal" = 18), name = NULL) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4 = ", x),
                                   iv2 = function(x) paste0("Interval = ", x))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = title, x = "b-parameter Mean (IV3)", y = "Reversal Rate (stacked by pair)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(size = 12))
}

# ── 분포별 그래프 (output/analysis/) ─────────────────────────────────────────
for (i in seq_along(IV4_LEVELS)) {
  dist <- IV4_LEVELS[i]
  df_d <- long_df %>% filter(iv4_theta_dist == dist)
  if (nrow(df_d) == 0) { message("건너뜀: ", dist, " 결과 없음"); next }

  title <- paste0("Ability Distribution: ", IV4_LABELS[i])
  ggsave(sprintf("output/analysis/transposition_pairs_bar_%s.png", dist),
         plot_bar(df_d, title), width = 11, height = 9, dpi = 150)
  ggsave(sprintf("output/analysis/transposition_pairs_line_%s.png", dist),
         plot_line(df_d, title), width = 11, height = 9, dpi = 150)
  ggsave(sprintf("output/analysis/transposition_pairs_stack_%s.png", dist),
         plot_stack(df_d, title), width = 11, height = 9, dpi = 150)
  message(sprintf("저장 완료: output/analysis/transposition_pairs_{bar,line,stack}_%s.png", dist))
}

# ── 정규분포 조건 그래프 (output/analysis_normal/) ───────────────────────────
normal_df <- long_df %>% filter(iv4_theta_dist == "normal")
if (nrow(normal_df) > 0) {
  dir.create("output/analysis_normal", recursive = TRUE, showWarnings = FALSE)
  ggsave("output/analysis_normal/transposition_pairs_normal_bar.png",
         plot_bar(normal_df, NULL), width = 11, height = 9, dpi = 150)
  ggsave("output/analysis_normal/transposition_pairs_normal.png",
         plot_line(normal_df, NULL), width = 11, height = 9, dpi = 150)
  ggsave("output/analysis_normal/transposition_pairs_normal_stack.png",
         plot_stack(normal_df, NULL), width = 11, height = 9, dpi = 150)
  message("저장 완료: output/analysis_normal/transposition_pairs_normal{,_bar,_stack}.png")
}
