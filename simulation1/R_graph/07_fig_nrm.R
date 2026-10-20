# =============================================================================
# 07_fig_nrm.R  (simulation1/R_graph)
# 결과 — NRM 경계모수 전치의 원인 (시뮬레이션2 결과를 임시로 사용)
#
# NRM 경계모수 b_k = (d_{k−1} − d_k) / (a1 · (ak_k − ak_{k−1}))
#   분모 = 추정 변별도 × 인접 채점함수 차이
#
#   fig14_nrm_denominator     그림 14  분모 a1·(ak4 − ak3) vs b̂4 (FIG14_COND 조건, 반복마다 점 1개)
#                                      빨강: 채점함수 전치 (ak4 ≤ ak3), 회색 점선: 분자를 참값으로 고정한 b4
#   fig15_nrm_cause           그림 15  전치 비율을 원인별로 나눈 누적 막대 (그림 8과 같은 배치·축)
#                                      극단 오차: b1~b4 중 하나라도 |b̂ − b| ≥ ERR_CUT
#                                      채점함수 이상: 인접 채점함수 차이 ≤ 0 또는 |차이| < DAK_CUT
#   extra_nrm_sf_reversal_crc 예비     채점함수가 전치된 반복 1개의 범주반응곡선 (참 모형 vs NRM 추정)
#
# 표
#   tab_nrm_cause_by_condition  조건별 원인별 전치 비율
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[07] NRM 전치 원인")

if (!file.exists(file.path(DATA_DIR, "sim2_nrm_item1.rds"))) {
  message("  시뮬레이션2 자료 없음 → 건너뜀")
} else {

reps2 <- prepare_reps(load_sim2()) %>% flag_nrm()
reps1 <- prepare_reps(load_sim1())

# =============================================================================
# 그림 14: 분모 vs b̂4
# =============================================================================
c14 <- reps2[match_cond(reps2, FIG14_COND), ] %>%
  mutate(denom = a1 * (ak4 - ak3),
         state = factor(ifelse(ak4 - ak3 <= 0, "rev", "ok"), levels = c("ok", "rev")))
if (nrow(c14) == 0) {
  message("  그림 14: 조건 자료 없음 → 건너뜀")
} else {
  sf4    <- FIG14_COND[["sf4"]]
  b4t    <- c14$b4_true[1]
  numer  <- DISCRIM * (4 - sf4) * b4t          # 참 분자 d3 − d4 = a·w4·b4
  out_n  <- sum(c14$b4 < FIG14_YLIM[1] | c14$b4 > FIG14_YLIM[2])
  message(sprintf("  그림 14: 반복 %d개 중 y축 범위(%g ~ %g) 밖 %d개, 채점함수 전치 %d개",
                  nrow(c14), FIG14_YLIM[1], FIG14_YLIM[2], out_n, sum(c14$state == "rev")))

  xr <- range(c(c14$denom, 0)); xr <- xr + c(-0.05, 0.05) * diff(xr)
  curve <- bind_rows(
    data.frame(x = seq(xr[1], -0.005, length.out = 300), part = "neg"),
    data.frame(x = seq(0.005, xr[2], length.out = 300), part = "pos")) %>%
    mutate(y = numer / x) %>% filter(y > FIG14_YLIM[1] - 5, y < FIG14_YLIM[2] + 5)

  p14 <- ggplot(c14, aes(denom, b4)) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.8) +
    geom_hline(yintercept = b4t, colour = "black", linetype = "dashed", linewidth = 0.9)
  if (isTRUE(FIG14_SHOW_CURVE))
    p14 <- p14 + geom_line(data = curve, aes(x, y, group = part), colour = "grey55",
                           linetype = "dotted", linewidth = LW)
  p14 <- p14 +
    geom_point(aes(colour = state), size = PS - 1, alpha = 0.7) +
    annotate("text", x = xr[2], y = b4t, label = "b[4]", parse = TRUE, vjust = -0.4, hjust = 1,
             size = BASE_SIZE / 3, family = FONT_FAMILY) +
    scale_colour_manual(values = c(ok = "grey35", rev = "#D7263D"),
                        labels = c(ok = "순서 유지", rev = "전치"),
                        name = LAB$sf_order, drop = FALSE) +
    coord_cartesian(xlim = xr, ylim = FIG14_YLIM) +
    labs(x = expression(hat(a)*(hat(s)[5] - hat(s)[4])), y = expression(hat(b)[4])) +
    theme_thesis() +
    guides(colour = guide_legend(override.aes = list(size = PS + 1, alpha = 1)))
  save_fig(p14, "fig14_nrm_denominator", 12, 9)
}

# =============================================================================
# 그림 15: 원인별 전치 비율 (그림 8과 같은 배치, 같은 y 범위)
# =============================================================================
cause_levels <- c("both", "extreme", "sf", "none")
cause_labels <- c(both = "극단 오차 + 채점함수 이상", extreme = "극단 오차만",
                  sf = "채점함수 이상만", none = "둘 다 아님")
cause_cols   <- c(both = "#7B3294", extreme = "#D55E00", sf = "#E69F00", none = "#999999")

cause <- reps2 %>%
  mutate(cause = case_when(extreme & sf_issue ~ "both", extreme ~ "extreme",
                           sf_issue ~ "sf", TRUE ~ "none")) %>%
  group_by(f_sf, f_int, f_mean) %>%
  summarise(n = n(),
            both    = sum(transposed & cause == "both") / n,
            extreme = sum(transposed & cause == "extreme") / n,
            sf      = sum(transposed & cause == "sf") / n,
            none    = sum(transposed & cause == "none") / n,
            .groups = "drop")
cause_long <- cause %>%
  pivot_longer(all_of(cause_levels), names_to = "cause", values_to = "prop") %>%
  mutate(cause = factor(cause, levels = cause_levels))

# 그림 8 (GPCM 쌍별), 예비 NRM 쌍별 그래프와 같은 y 범위
s1 <- summarise_transposition(reps1); s2 <- summarise_transposition(reps2)
stack_top <- nice_upper(unlist(lapply(list(s1, s2), function(s)
  pmax(s$prop_any, s$prop_12 + s$prop_23 + s$prop_34))))

p15 <- plot_stack(cause_long, "cause", cause_cols, cause_labels, LAB$cause, stack_top) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))
save_fig(p15, "fig15_nrm_cause", 15, 13)

save_table(cause %>% transmute(`채점함수 간격` = f_sf, `경계모수 간격` = f_int, `문항 난이도` = f_mean,
                               `반복 수` = n,
                               `전체 전치(%)` = round(100 * (both + extreme + sf + none), 1),
                               `극단 오차 + 채점함수 이상(%)` = round(100 * both, 1),
                               `극단 오차만(%)` = round(100 * extreme, 1),
                               `채점함수 이상만(%)` = round(100 * sf, 1),
                               `둘 다 아님(%)` = round(100 * none, 1)),
           "tab_nrm_cause_by_condition")

# =============================================================================
# 예비: 채점함수가 전치된 반복 1개의 범주반응곡선 (FIG14_COND 조건)
#   ak4 ≤ ak3 인 반복 중 (ak4 − ak3) 가 그 반복들의 중앙값에 가장 가까운 반복
# =============================================================================
rev_reps <- reps2[match_cond(reps2, FIG14_COND) & reps2$sf_rev4, ]
if (nrow(rev_reps) == 0) {
  message("  예비 범주반응곡선: 채점함수가 전치된 반복이 없음 → 건너뜀")
} else {
  gap <- rev_reps$ak4 - rev_reps$ak3
  r   <- rev_reps[which.min(abs(gap - median(gap))), ]
  th  <- seq(-4, 6, length.out = 400)
  b_t <- unlist(r[paste0("b", 1:4, "_true")])
  panel <- c(LAB$true_model, LAB$nrm)
  dd <- bind_rows(
    crc_long(th, crc_true(th, b_t, sf = c(0, 1, 2, FIG14_COND[["sf4"]], 4)), panel = panel[1]),
    crc_long(th, crc_nrm(th, r$a1, unlist(r[paste0("ak", 0:4)]), unlist(r[paste0("d", 0:4)])),
             panel = panel[2]))
  dd$panel <- factor(dd$panel, levels = panel)
  lab <- crc_peak_labels(dd, "panel")
  p_ex <- ggplot(dd, aes(theta, prob, colour = category)) +
    geom_line(linewidth = LW) +
    geom_label(data = lab, aes(label = category), vjust = -0.25, size = BASE_SIZE / 3.2,
               family = FONT_FAMILY, label.size = 0, fill = alpha("white", 0.8), show.legend = FALSE) +
    scale_colour_manual(values = COL_CATEGORY, guide = "none") +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
    facet_wrap(~panel, nrow = 1) +
    labs(x = LAB$theta, y = LAB$prob) +
    theme_thesis()
  save_fig(p_ex, "extra_nrm_sf_reversal_crc", 15, 6.5)
  message(sprintf("  예비 범주반응곡선: 반복 %d (ak3 = %.3f, ak4 = %.3f, b4 = %.3f)",
                  r$rep_id, r$ak3, r$ak4, r$b4))
}

}
