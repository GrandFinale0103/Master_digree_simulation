# =============================================================================
# 14_plot_nrm_transposition_pairs.R  (simulation2, NRM 추정 결과)
# NRM 경계모수 전치가 어느 쌍(b1&b2, b2&b3, b3&b4)에서 발생했는지 시각화
#
# 배치: 원래 전치 비율 그래프와 동일
#   행 = IV1(s4), 열 = IV2(간격, 왼쪽 좁음 → 오른쪽 넓음), 패널 x축 = IV3(평균)
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력)
#   전치 판정: b_k ≥ b_{k+1}  (수렴한 반복만, 비율 분모 = 수렴 반복 수)
#   한 반복에서 여러 쌍이 동시에 전치될 수 있음
#
# 출력 (output/nrm/analysis/):
#   nrm_transposition_pairs_summary.csv — 조건별 쌍별 전치 수·비율
#   nrm_transposition_pairs_bar.png     — 쌍별 막대 그래프
#   nrm_transposition_pairs_line.png    — 쌍별 선 그래프
#   nrm_transposition_pairs_stack.png   — 쌍별 누적 막대 + 전체 전치 비율(검은 점)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

OUT_DIR <- "output/nrm/analysis"
in_path <- file.path(OUT_DIR, "nrm_test_rep.csv")
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)

make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

PAIR_LEVELS  <- c("prop_12", "prop_23", "prop_34")
PAIR_LABELS  <- c("b1 & b2", "b2 & b3", "b3 & b4")
pair_colours <- c("b1 & b2" = "#4C9BE8", "b2 & b3" = "#F2A541", "b3 & b4" = "#D1495B")

rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character"))

pairs_df <- rep_df %>%
  filter(as.logical(converged) %in% TRUE, !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4)) %>%
  mutate(t12 = b1 >= b2, t23 = b2 >= b3, t34 = b3 >= b4,
         t_any = t12 | t23 | t34) %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
  summarise(n_converged  = n(),
            n_transposed = sum(t_any),
            n_12 = sum(t12), n_23 = sum(t23), n_34 = sum(t34),
            .groups = "drop") %>%
  mutate(prop_transposed = n_transposed / n_converged,
         prop_12 = n_12 / n_converged,
         prop_23 = n_23 / n_converged,
         prop_34 = n_34 / n_converged)

write.csv(pairs_df, file.path(OUT_DIR, "nrm_transposition_pairs_summary.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s/nrm_transposition_pairs_summary.csv (%d 조건)",
                OUT_DIR, nrow(pairs_df)))

long_df <- pairs_df %>%
  pivot_longer(all_of(PAIR_LEVELS), names_to = "pair", values_to = "prop") %>%
  mutate(
    pair = factor(pair, levels = PAIR_LEVELS, labels = PAIR_LABELS),
    iv1  = make_num_factor(iv1_sf4,        "%.2f"),
    iv2  = make_num_factor(iv2_b_interval, "%.2f", decreasing = FALSE),  # 좁은 → 넓은
    iv3  = make_num_factor(iv3_b_mean,     "%.1f")
  )

common_layers <- list(
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                 iv2 = function(x) paste0("Int=", x))),
  labs(x = "b-parameter Mean (IV3)"),
  theme_bw(base_size = 9),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
pct_0_1 <- scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))

save_plot <- function(p, fname) {
  ggsave(file.path(OUT_DIR, fname), p, width = 24, height = 18, dpi = 150)
  message(sprintf("저장 완료: %s/%s", OUT_DIR, fname))
}

# ── 막대 ─────────────────────────────────────────────────────────────────────
p_bar <- ggplot(long_df, aes(x = iv3, y = prop, fill = pair)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = pair_colours, name = "Reversed pair") +
  pct_0_1 + labs(y = "Reversal Rate (NRM)") + common_layers
save_plot(p_bar, "nrm_transposition_pairs_bar.png")

# ── 선 ───────────────────────────────────────────────────────────────────────
p_line <- ggplot(long_df, aes(x = iv3, y = prop, colour = pair, shape = pair, group = pair)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  scale_colour_manual(values = pair_colours, name = "Reversed pair") +
  scale_shape_discrete(name = "Reversed pair") +
  pct_0_1 + labs(y = "Reversal Rate (NRM)") + common_layers
save_plot(p_line, "nrm_transposition_pairs_line.png")

# ── 누적 막대 + 전체 전치 비율(점) ───────────────────────────────────────────
total_df <- long_df %>% distinct(iv1, iv2, iv3, prop_transposed)
p_stack <- ggplot(long_df, aes(x = iv3, y = prop)) +
  geom_col(aes(fill = pair), position = position_stack(reverse = TRUE), width = 0.7) +
  geom_point(data = total_df, aes(y = prop_transposed, shape = "Any reversal"),
             size = 1.6, colour = "black") +
  scale_fill_manual(values = pair_colours, name = "Reversed pair") +
  scale_shape_manual(values = c("Any reversal" = 18), name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(y = "Reversal Rate (NRM, stacked by pair)") + common_layers
save_plot(p_stack, "nrm_transposition_pairs_stack.png")
