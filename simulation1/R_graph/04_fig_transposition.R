# =============================================================================
# 04_fig_transposition.R  (simulation1/R_graph)
# 결과 — 경계모수 전치 비율 (시뮬레이션1 GPCM, 시뮬레이션2 NRM)
#
# 전치 = 인접 경계모수 쌍 중 하나라도 b_k ≥ b_{k+1}
# 분모 = 수렴하고 b1~b4 가 모두 계산된 반복
#
# 배치 (모두 같음): 9개 패널 막대 그래프
#   행 = 채점함수 간격, 열 = 경계모수 간격, x = 문항 난이도
# y축: 그림 7·8·12·15 와 예비 NRM 쌍별 그림이 모두 같음 (trans_ylim)
#
#   fig07_sim1_transposition         그림 7   GPCM 전치 비율
#   fig12_sim2_transposition         그림 12  NRM 전치 비율
#   fig12b_transposition_difference  그림 12  NRM − GPCM 전치 비율 차이
#   fig08_sim1_transposition_pairs   그림 8   GPCM 쌍별 전치 비율 누적 막대
#   extra_sim2_transposition_pairs   예비     NRM 쌍별 전치 비율 누적 막대
#
# 표
#   tab_transposition_by_condition   조건별 전치 비율 (두 시뮬레이션)
#   tab_sim1_transposition_lrt       GPCM 전치 여부 로지스틱 회귀 우도비 검정
#   tab_sim2_transposition_lrt       NRM 전치 여부 로지스틱 회귀 우도비 검정
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[04] 전치 비율")

reps1 <- prepare_reps(load_sim1())
sum1  <- summarise_transposition(reps1)

has2 <- file.exists(file.path(DATA_DIR, "sim2_nrm_item1.rds"))
if (has2) {
  reps2 <- prepare_reps(load_sim2())
  sum2  <- summarise_transposition(reps2)
} else {
  message("  시뮬레이션2 자료 없음 → NRM 그래프 건너뜀")
}
ylim_t <- trans_ylim(if (has2) list(sum1, sum2) else list(sum1))

W9 <- 15; H9 <- 13   # 9개 패널 그림 크기 (인치)

# ── 그림 7·12: 전치 비율 ────────────────────────────────────────────────────
save_fig(plot_bars(sum1, "prop_any", LAB$trans, ylim_t), "fig07_sim1_transposition", W9, H9)
if (has2) {
  save_fig(plot_bars(sum2, "prop_any", LAB$trans, ylim_t), "fig12_sim2_transposition", W9, H9)

  diff_df <- sum1 %>% select(f_sf, f_int, f_mean, gpcm = prop_any) %>%
    inner_join(sum2 %>% select(f_sf, f_int, f_mean, nrm = prop_any), by = c("f_sf", "f_int", "f_mean")) %>%
    mutate(diff = nrm - gpcm)
  lim <- c(-1, 1) * nice_upper(abs(diff_df$diff))
  save_fig(plot_bars(diff_df, "diff", LAB$trans_diff, lim, fill_by_sign = TRUE,
                     y_labels = scales::label_number(accuracy = 1, scale = 100, suffix = "%p")),
           "fig12b_transposition_difference", W9, H9)
}

# ── 그림 8: 쌍별 전치 비율 누적 막대 ────────────────────────────────────────
pair_long <- function(s) {
  s %>% pivot_longer(c(prop_12, prop_23, prop_34), names_to = "pair",
                     names_prefix = "prop_", values_to = "prop") %>%
    mutate(pair = factor(pair, levels = c("12", "23", "34")))
}
pair_labels <- c(expression(group("(", list(b[1], b[2]), ")")),
                 expression(group("(", list(b[2], b[3]), ")")),
                 expression(group("(", list(b[3], b[4]), ")")))
plot_pairs <- function(s) {
  plot_stack(pair_long(s), "pair", COL_PAIR, pair_labels, LAB$pair, ylim_t,
             total = s, total_label = LAB$trans_any)
}
save_fig(plot_pairs(sum1), "fig08_sim1_transposition_pairs", W9, H9 + 1)
if (has2) save_fig(plot_pairs(sum2), "extra_sim2_transposition_pairs", W9, H9 + 1)

# ── 표 ──────────────────────────────────────────────────────────────────────
tab_cond <- function(s, model) {
  s %>% transmute(`모형` = model, `채점함수 간격` = f_sf, `경계모수 간격` = f_int,
                  `문항 난이도` = f_mean, `반복 수` = n,
                  `전치 비율(%)` = round(100 * prop_any, 1),
                  `b1·b2 전치(%)` = round(100 * prop_12, 1),
                  `b2·b3 전치(%)` = round(100 * prop_23, 1),
                  `b3·b4 전치(%)` = round(100 * prop_34, 1))
}
save_table(bind_rows(tab_cond(sum1, "GPCM"), if (has2) tab_cond(sum2, "NRM")),
           "tab_transposition_by_condition")
save_table(transposition_lrt(reps1), "tab_sim1_transposition_lrt")
if (has2) save_table(transposition_lrt(reps2), "tab_sim2_transposition_lrt")
