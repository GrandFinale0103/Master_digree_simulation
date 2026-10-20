# =============================================================================
# 03_fig_design.R  (simulation1/R_graph)
# 방법 — 설계 그림 (00_config.R 의 시뮬레이션1 조건으로 그림)
#
#   fig06_true_boundaries  그림 6  조건별 참 경계모수 위치와 θ ~ N(0, 1) 분포
#                                  행 = 경계모수 간격, 열 = 문항 심각도
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[03] 설계 그림")

theta <- seq(-4, 4.5, length.out = 400)
grid  <- expand.grid(interval = SIM1_IV2, mean = SIM1_IV3, sf4 = SIM1_IV1[1])

dens <- grid %>%
  tidyr::crossing(theta = theta) %>%
  mutate(dens = dnorm(theta)) %>%
  add_design_factors()

pts <- grid %>%
  tidyr::crossing(k = 1:4) %>%
  mutate(b = mean + interval * B_OFFSET[k], lab = paste0("b[", k, "]"),
         h = ifelse(k %% 2 == 1, 0.92, 0.78)) %>%   # 라벨이 겹치지 않게 높이를 번갈아
  add_design_factors()

top <- max(dnorm(0)) * 1.25
p6 <- ggplot() +
  geom_area(data = dens, aes(theta, dens), fill = "grey85", colour = "grey55", linewidth = 0.6) +
  geom_segment(data = pts, aes(x = b, xend = b, y = 0, yend = top * (h - 0.06), colour = factor(k)),
               linewidth = LW) +
  geom_point(data = pts, aes(b, 0, colour = factor(k)), size = PS) +
  geom_text(data = pts, aes(b, top * h, label = lab, colour = factor(k)), parse = TRUE,
            size = BASE_SIZE / 3.4, family = FONT_FAMILY) +
  scale_colour_manual(values = unname(COL_PARAM), guide = "none") +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(limits = c(0, top), breaks = NULL, expand = expansion(mult = c(0, 0.02))) +
  facet_grid(rows = vars(f_int), cols = vars(f_mean),
             labeller = labeller(f_int = lab_interval, f_mean = lab_mean)) +
  labs(x = LAB$theta, y = NULL) +
  theme_thesis()
save_fig(p6, "fig06_true_boundaries", 16, 10)
