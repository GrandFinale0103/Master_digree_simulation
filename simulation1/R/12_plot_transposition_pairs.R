# =============================================================================
# 12_plot_transposition_pairs.R  (simulation1, GPCM 추정 결과)
# 전치가 어느 경계모수 쌍(b1&b2, b2&b3, b3&b4)에서 발생했는지 시각화
#
# 입력: output/analysis/transposition_summary.csv  (08_analysis_transposition.R 출력)
#   n_12 / n_23 / n_34 : 해당 쌍이 전치된 반복 수 (한 반복에서 여러 쌍 중복 가능)
#   비율 분모          : 수렴한 반복 수 (n_converged)
#
# 출력:
#   output/analysis/transposition_pairs_plot.png         — 행 IV4×IV1, 열 IV2, x축 IV3
#   output/analysis/transposition_pairs_by_interval.png  — x축 IV2(간격), 행 IV3, 열 IV4×IV1
#   output/analysis/transposition_pairs_summary.csv      — 그래프에 쓴 비율 표
#   output/analysis_normal/transposition_pairs_normal.png — 정규분포 조건만 (선 그래프)
#                                                          x축 IV2, 행 IV3, 열 IV1
#   output/analysis_normal/transposition_pairs_normal_bar.png — 같은 내용의 막대 그래프
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

PAIR_LEVELS <- c("prop_12", "prop_23", "prop_34")
PAIR_LABELS <- c("b1 & b2", "b2 & b3", "b3 & b4")
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
    iv3  = make_num_factor(iv3_b_mean,     "%.1f"),
    iv4  = factor(iv4_theta_dist, levels = IV4_LEVELS, labels = IV4_LABELS)
  )

# ── 그래프 1: 조건별 쌍별 전치 비율 막대 ─────────────────────────────────────
p1 <- ggplot(long_df, aes(x = iv3, y = prop, fill = pair)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_grid(rows = vars(iv4, iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4 = ", x),
                                 iv2 = function(x) paste0("Interval = ", x),
                                 iv4 = label_value)) +
  scale_fill_manual(values = pair_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "b-parameter Mean (IV3)", y = "Reversal Rate", fill = "Reversed pair") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text.y = element_text(size = 8, angle = 0),
        strip.text.x = element_text(size = 10))

ggsave("output/analysis/transposition_pairs_plot.png", p1,
       width = 12, height = 20, dpi = 150)
message("저장 완료: output/analysis/transposition_pairs_plot.png")

# ── 그래프 2: x축 = 경계모수 간격, 쌍별 선 ────────────────────────────────────
p2 <- ggplot(long_df, aes(x = iv2_b_interval, y = prop, colour = pair)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_grid(rows = vars(iv3), cols = vars(iv4, iv1),
             labeller = labeller(iv3 = function(x) paste0("Mean = ", x),
                                 iv1 = function(x) paste0("s4 = ", x),
                                 iv4 = label_value)) +
  scale_colour_manual(values = pair_colours) +
  scale_x_continuous(breaks = sort(unique(long_df$iv2_b_interval))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "b-parameter Interval (IV2)", y = "Reversal Rate", colour = "Reversed pair") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

ggsave("output/analysis/transposition_pairs_by_interval.png", p2,
       width = 24, height = 8, dpi = 150)
message("저장 완료: output/analysis/transposition_pairs_by_interval.png")

# ── 그래프 3: 정규분포 조건만 (27 조건) ─────────────────────────────────────
# 간격 수준이 3개뿐이라 막대보다 선이 간격에 따른 감소 추세를 잘 보여줌
dir.create("output/analysis_normal", recursive = TRUE, showWarnings = FALSE)
normal_df <- long_df %>% filter(iv4_theta_dist == "normal")

if (nrow(normal_df) == 0) {
  message("건너뜀: 정규분포 조건 결과가 없습니다.")
} else {
  p3 <- ggplot(normal_df, aes(x = iv2_b_interval, y = prop,
                              colour = pair, shape = pair)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 3) +
    facet_grid(rows = vars(iv3), cols = vars(iv1),
               labeller = labeller(iv3 = function(x) paste0("b Mean = ", x),
                                   iv1 = function(x) paste0("s4 = ", x))) +
    scale_colour_manual(values = pair_colours) +
    scale_x_continuous(breaks = sort(unique(normal_df$iv2_b_interval))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = "b-parameter Interval (IV2)", y = "Reversal Rate",
         colour = "Reversed pair", shape = "Reversed pair") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(size = 12))

  ggsave("output/analysis_normal/transposition_pairs_normal.png", p3,
         width = 11, height = 9, dpi = 150)
  message("저장 완료: output/analysis_normal/transposition_pairs_normal.png")
}

# ── 그래프 4: 정규분포 조건만 — 막대 그래프 ─────────────────────────────────
if (nrow(normal_df) > 0) {
  p4 <- ggplot(normal_df, aes(x = iv2, y = prop, fill = pair)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    facet_grid(rows = vars(iv3), cols = vars(iv1),
               labeller = labeller(iv3 = function(x) paste0("b Mean = ", x),
                                   iv1 = function(x) paste0("s4 = ", x))) +
    scale_fill_manual(values = pair_colours) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = "b-parameter Interval (IV2)", y = "Reversal Rate", fill = "Reversed pair") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(size = 12))

  ggsave("output/analysis_normal/transposition_pairs_normal_bar.png", p4,
         width = 11, height = 9, dpi = 150)
  message("저장 완료: output/analysis_normal/transposition_pairs_normal_bar.png")
}
