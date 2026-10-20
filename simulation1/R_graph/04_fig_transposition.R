# =============================================================================
# 04_fig_transposition.R  (simulation1/R_graph)
# 결과 — 경계모수 전치 비율 (시뮬레이션1 GPCM, 시뮬레이션2 NRM)
#
# 전치 = 인접 경계모수 쌍 중 하나라도 b_k ≥ b_{k+1}
# 분모 = 수렴하고 b1~b4 가 모두 계산된 반복
#
# 비교되는 그래프는 축·패널 배열이 같음
#   fig07_sim1_transposition         그림 7   GPCM 전치 비율 (열 = 채점함수 간격, x = 난이도, 색 = 경계모수 간격)
#   fig12_sim2_transposition         그림 12  NRM 전치 비율  (그림 7과 같은 축)
#   fig12b_transposition_difference  그림 12  NRM − GPCM 전치 비율 차이 (같은 배치)
#   fig08_sim1_transposition_pairs   그림 8   GPCM 쌍별 전치 비율 누적 막대 (행 = 채점함수, 열 = 경계모수 간격)
#   extra_sim2_transposition_pairs   예비     NRM 쌍별 전치 비율 (그림 8과 같은 축)
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
all_sum <- if (has2) list(sum1, sum2) else list(sum1)

# ── 그림 7·12: 전치 비율 (같은 y 범위) ──────────────────────────────────────
y_top <- nice_upper(unlist(lapply(all_sum, `[[`, "prop_any")))
pct   <- scales::label_percent(accuracy = 1)

save_fig(plot_lines(sum1, "prop_any", LAB$trans, ylim = c(0, y_top), y_labels = pct),
         "fig07_sim1_transposition", 16, 7)
if (has2) {
  save_fig(plot_lines(sum2, "prop_any", LAB$trans, ylim = c(0, y_top), y_labels = pct),
           "fig12_sim2_transposition", 16, 7)

  diff_df <- sum1 %>% select(f_sf, f_int, f_mean, gpcm = prop_any) %>%
    inner_join(sum2 %>% select(f_sf, f_int, f_mean, nrm = prop_any), by = c("f_sf", "f_int", "f_mean")) %>%
    mutate(diff = nrm - gpcm)
  lim <- max(abs(diff_df$diff)); lim <- c(-1, 1) * nice_upper(lim)
  save_fig(plot_lines(diff_df, "diff", LAB$trans_diff, ylim = lim, hline = 0,
                      y_labels = scales::label_number(accuracy = 1, scale = 100, suffix = "%p")),
           "fig12b_transposition_difference", 16, 7)
}

# ── 그림 8: 쌍별 전치 비율 누적 막대 (같은 y 범위) ──────────────────────────
pair_long <- function(s) {
  s %>% pivot_longer(c(prop_12, prop_23, prop_34), names_to = "pair",
                     names_prefix = "prop_", values_to = "prop") %>%
    mutate(pair = factor(pair, levels = c("12", "23", "34")))
}
pair_labels <- c(expression(group("(", list(b[1], b[2]), ")")),
                 expression(group("(", list(b[2], b[3]), ")")),
                 expression(group("(", list(b[3], b[4]), ")")))
stack_top <- nice_upper(unlist(lapply(all_sum, function(s)
  pmax(s$prop_any, s$prop_12 + s$prop_23 + s$prop_34))))

plot_pairs <- function(s) {
  plot_stack(pair_long(s), "pair", COL_PAIR, pair_labels, LAB$pair, stack_top,
             total = s, total_label = LAB$trans_any)
}
save_fig(plot_pairs(sum1), "fig08_sim1_transposition_pairs", 15, 13)
if (has2) save_fig(plot_pairs(sum2), "extra_sim2_transposition_pairs", 15, 13)

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
