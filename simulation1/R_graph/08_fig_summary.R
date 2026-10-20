# =============================================================================
# 08_fig_summary.R  (simulation1/R_graph)
# 논의 — 27개 조건의 주요 결과를 한 장에 요약한 히트맵
#
#   fig16_summary_heatmap   행 = 조건 (채점함수 간격별로 묶음), 열 = 지표
#     GPCM 전치 비율 · GPCM b4 편향 · NRM 전치 비율 · NRM b4 RMSE(이상치 제거 후)
#     칸 색 = 열 안에서 |값| / 최댓값 (진할수록 큼), 칸 숫자 = 실제 값
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[08] 요약 히트맵")

reps1 <- prepare_reps(load_sim1())
t1 <- summarise_transposition(reps1)
e1 <- summarise_error(reps1) %>% filter(param == "b4")

cols <- list(
  t1 %>% transmute(f_sf, f_int, f_mean, key = "gpcm_trans", value = prop_any),
  e1 %>% transmute(f_sf, f_int, f_mean, key = "gpcm_bias", value = bias))
col_lab <- c(gpcm_trans = "GPCM\n전치 비율", gpcm_bias = "GPCM\nb4 편향",
             nrm_trans  = "NRM\n전치 비율",  nrm_rmse  = "NRM\nb4 RMSE")

if (file.exists(file.path(DATA_DIR, "sim2_nrm_item1.rds"))) {
  reps2 <- prepare_reps(load_sim2())
  t2 <- summarise_transposition(reps2)
  e2 <- summarise_error(reps2, exclude = nrm_outlier_matrix(reps2)) %>% filter(param == "b4")
  cols <- c(cols, list(
    t2 %>% transmute(f_sf, f_int, f_mean, key = "nrm_trans", value = prop_any),
    e2 %>% transmute(f_sf, f_int, f_mean, key = "nrm_rmse", value = rmse)))
}

hm <- bind_rows(cols) %>%
  mutate(key = factor(key, levels = intersect(names(col_lab), unique(key))),
         row = factor(paste0(LAB$interval, " ", f_int, " · ", LAB$mean, " ", f_mean),
                      levels = rev(as.vector(t(outer(paste0(LAB$interval, " ", INTERVAL_LEVELS),
                                                     paste0(LAB$mean, " ", MEAN_LEVELS),
                                                     function(a, b) paste0(a, " · ", b))))))) %>%
  group_by(key) %>% mutate(intensity = abs(value) / max(abs(value), na.rm = TRUE)) %>% ungroup() %>%
  mutate(label = ifelse(key %in% c("gpcm_trans", "nrm_trans"),
                        sprintf("%.0f%%", 100 * value), sub("^-(0\\.0+)$", "\\1", sprintf("%.2f", value))))

p16 <- ggplot(hm, aes(key, row, fill = intensity)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = label, colour = intensity > 0.6), size = BASE_SIZE / 3.6,
            family = FONT_FAMILY, show.legend = FALSE) +
  scale_fill_gradient(low = "#FFF5EB", high = "#A63603", guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_x_discrete(labels = col_lab, position = "top") +
  facet_grid(rows = vars(f_sf), scales = "free_y", space = "free_y",
             labeller = labeller(f_sf = lab_sf)) +
  labs(x = NULL, y = NULL) +
  theme_thesis() +
  theme(panel.grid = element_blank(), strip.text.y = element_text(angle = 0),
        axis.ticks = element_blank())
save_fig(p16, "fig16_summary_heatmap", 16, 15)
