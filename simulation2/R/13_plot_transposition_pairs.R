# =============================================================================
# 13_plot_transposition_pairs.R  (simulation2, GPCM 추정 결과)
# 전치가 어느 경계모수 쌍(b1&b2, b2&b3, b3&b4)에서 발생했는지 시각화
#
# 목적: 경계모수 간격이 좁을 때의 전치가 특정 쌍에 몰리지 않고 세 쌍 모두에서
#       비슷하게 나타나는지(= 추정 오차), 간격이 넓어져도 b3&b4 에만 남는지
#       (= 채점함수 등 체계적 원인) 비교
#
# 입력: output/analysis/transposition_summary.csv  (08_analysis_transposition.R 출력)
#   n_12 / n_23 / n_34 : 해당 쌍이 전치된 반복 수 (한 반복에서 여러 쌍 중복 가능)
#   비율 분모          : 수렴한 반복 수 (n_converged)
#
# 출력:
#   output/analysis/transposition_pairs_plot.png         — 조건별 쌍별 전치 비율
#                                                          (행 IV1 × 열 IV2, x축 IV3)
#   output/analysis/transposition_pairs_line.png         — 같은 배치의 쌍별 선 그래프
#   output/analysis/transposition_pairs_summary.csv      — 그래프에 쓴 비율 표
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

PAIR_LEVELS <- c("prop_12", "prop_23", "prop_34")
PAIR_LABELS <- c("b1 & b2", "b2 & b3", "b3 & b4")

pairs_df <- read.csv(in_path, stringsAsFactors = FALSE) %>%
  mutate(
    prop_12 = n_12 / n_converged,
    prop_23 = n_23 / n_converged,
    prop_34 = n_34 / n_converged
  )

write.csv(pairs_df %>%
            select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, n_converged,
                   n_transposed, prop_transposed, n_12, n_23, n_34,
                   prop_12, prop_23, prop_34),
          "output/analysis/transposition_pairs_summary.csv", row.names = FALSE)

long_df <- pairs_df %>%
  pivot_longer(all_of(PAIR_LEVELS), names_to = "pair", values_to = "prop") %>%
  mutate(
    pair = factor(pair, levels = PAIR_LEVELS, labels = PAIR_LABELS),
    iv1  = make_num_factor(iv1_sf4,        "%.2f"),
    iv2  = make_num_factor(iv2_b_interval, "%.2f", decreasing = FALSE),  # 좁은 → 넓은
    iv3  = make_num_factor(iv3_b_mean,     "%.1f")
  )

pair_colours <- c("b1 & b2" = "#4C9BE8", "b2 & b3" = "#F2A541", "b3 & b4" = "#D1495B")

# ── 그래프 1: 기존 전치 비율 그래프와 같은 배치, 쌍별 막대 ────────────────────
p1 <- ggplot(long_df, aes(x = iv3, y = prop, fill = pair)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                 iv2 = function(x) paste0("Int=", x))) +
  scale_fill_manual(values = pair_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "b-parameter Mean (IV3)", y = "Reversal Rate", fill = "Reversed pair") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/analysis/transposition_pairs_plot.png", p1,
       width = 24, height = 18, dpi = 150)
message("저장 완료: output/analysis/transposition_pairs_plot.png")

# ── 그래프 2: 같은 배치(행 IV1 × 열 IV2, x축 IV3), 쌍별 선 ─────────────────────
p2 <- ggplot(long_df, aes(x = iv3, y = prop, colour = pair, shape = pair, group = pair)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                 iv2 = function(x) paste0("Int=", x))) +
  scale_colour_manual(values = pair_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "b-parameter Mean (IV3)", y = "Reversal Rate",
       colour = "Reversed pair", shape = "Reversed pair") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("output/analysis/transposition_pairs_line.png", p2,
       width = 24, height = 18, dpi = 150)
message("저장 완료: output/analysis/transposition_pairs_line.png")
