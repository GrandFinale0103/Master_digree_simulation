# =============================================================================
# 05_fig_error.R  (simulation1/R_graph)
# 결과 — 경계모수 추정 편향·RMSE
#
# 편향 = mean(b̂ − b),  RMSE = sqrt(mean((b̂ − b)²))  (수렴하고 b1~b4 가 계산된 반복)
# 시뮬레이션2(NRM)는 이상치 제거 후 계산
#   이상치: |b̂ − b| ≥ ERR_CUT (00_config.R), OUTLIER_UNIT 에 따라 반복 전체 또는 해당 값만 제거
#
# 배치 (모든 그림 같음): 행 = 경계모수, 열 = 채점함수 간격, x = 문항 난이도, 색 = 경계모수 간격
# 같은 지표(편향끼리, RMSE끼리)는 시뮬레이션1·2의 행별 y 범위를 같게 맞춤
#   *_b3b4 : b3, b4 만 (발표용)        *_b1b4 : b1~b4 (부록용)
#
#   fig09_sim1_bias_*   그림 9   GPCM 편향
#   fig13_sim2_rmse_*   그림 13  NRM RMSE (이상치 제거 후)
#   extra_sim1_rmse_*   예비     GPCM RMSE      (그림 13과 같은 축)
#   extra_sim2_bias_*   예비     NRM 편향 (이상치 제거 후, 그림 9와 같은 축)
#
# 표
#   tab_sim1_bias_anova_b3 / _b4   GPCM (b̂ − b) Type III 분산분석
#   tab_error_by_condition         조건 × 경계모수별 편향·RMSE (두 시뮬레이션)
#   tab_sim2_outliers              NRM 조건별 이상치 제거 현황
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[05] 편향·RMSE")

PARAM_SETS <- list(b3b4 = c("b3", "b4"), b1b4 = c("b1", "b2", "b3", "b4"))
add_param_lab <- function(df) {
  df %>% mutate(param_lab = factor(paste0("b[", sub("b", "", param), "]"),
                                   levels = paste0("b[", 1:4, "]")))
}

reps1 <- prepare_reps(load_sim1())
err1  <- summarise_error(reps1) %>% add_param_lab()

has2 <- file.exists(file.path(DATA_DIR, "sim2_nrm_item1.rds"))
if (has2) {
  reps2 <- prepare_reps(load_sim2())
  out2  <- nrm_outlier_matrix(reps2)
  err2  <- summarise_error(reps2, exclude = out2) %>% add_param_lab()
  err2_all <- summarise_error(reps2) %>% add_param_lab()
}

# ── 그림 ────────────────────────────────────────────────────────────────────
height_for <- function(params) 4 + 3.2 * length(params)

for (set in names(PARAM_SETS)) {
  ps <- PARAM_SETS[[set]]
  e1 <- err1 %>% filter(param %in% ps) %>% droplevels()
  e2 <- if (has2) err2 %>% filter(param %in% ps) %>% droplevels()
  both <- if (has2) list(e1, e2) else list(e1)
  h <- height_for(ps)

  bias_blank <- row_range_blank(both, "bias", rows = "param_lab", include_zero = TRUE)
  rmse_blank <- row_range_blank(both, "rmse", rows = "param_lab", include_zero = TRUE)

  save_fig(plot_lines(e1, "bias", LAB$bias, rows = "param_lab", blank = bias_blank, hline = 0),
           paste0("fig09_sim1_bias_", set), 16, h)
  save_fig(plot_lines(e1, "rmse", LAB$rmse, rows = "param_lab", blank = rmse_blank),
           paste0("extra_sim1_rmse_", set), 16, h)
  if (has2) {
    save_fig(plot_lines(e2, "rmse", LAB$rmse, rows = "param_lab", blank = rmse_blank),
             paste0("fig13_sim2_rmse_", set), 16, h)
    save_fig(plot_lines(e2, "bias", LAB$bias, rows = "param_lab", blank = bias_blank, hline = 0),
             paste0("extra_sim2_bias_", set), 16, h)
  }
}

# ── 표: 분산분석 (GPCM b3, b4) ──────────────────────────────────────────────
save_table(anova_type3(reps1, "err_b3"), "tab_sim1_bias_anova_b3")
save_table(anova_type3(reps1, "err_b4"), "tab_sim1_bias_anova_b4")

# ── 표: 조건별 편향·RMSE ────────────────────────────────────────────────────
tab_err <- function(e, model) {
  e %>% transmute(`모형` = model, `경계모수` = param, `채점함수 간격` = f_sf,
                  `경계모수 간격` = f_int, `문항 난이도` = f_mean, `반복 수` = n,
                  `편향` = round(bias, 3), RMSE = round(rmse, 3))
}
save_table(bind_rows(tab_err(err1, "GPCM"),
                     if (has2) tab_err(err2, "NRM(이상치 제거 후)"),
                     if (has2) tab_err(err2_all, "NRM(전체)")),
           "tab_error_by_condition")

# ── 표: NRM 이상치 제거 현황 ────────────────────────────────────────────────
if (has2) {
  rep_out <- apply(out2, 1, any)
  out_tab <- reps2 %>% mutate(.out = rep_out) %>%
    group_by(f_sf, f_int, f_mean) %>%
    summarise(n = n(), n_out = sum(.out), .groups = "drop")
  out_tab <- bind_rows(
    out_tab %>% mutate(across(c(f_sf, f_int, f_mean), as.character)),
    data.frame(f_sf = "전체", f_int = "", f_mean = "", n = sum(out_tab$n), n_out = sum(out_tab$n_out))) %>%
    transmute(`채점함수 간격` = f_sf, `경계모수 간격` = f_int, `문항 난이도` = f_mean,
              `분석 반복 수` = n, `제거 반복 수` = n_out,
              `제거 비율(%)` = round(100 * n_out / n, 1))
  save_table(out_tab, "tab_sim2_outliers")
  message(sprintf("  NRM 이상치 기준 |오차| ≥ %g (%s 단위): 전체 %d / %d 반복에 이상치",
                  ERR_CUT, ifelse(OUTLIER_UNIT == "rep", "반복", "값"), sum(rep_out), length(rep_out)))
}
