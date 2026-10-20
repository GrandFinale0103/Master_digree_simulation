# =============================================================================
# 02_fig_concept.R  (simulation1/R_graph)
# 연구 목적 — 개념 그림 (자료 없이 00_config.R 의 예시 모수로 그림)
#
#   fig01_crc_discrimination   그림 1   변별도에 따른 범주반응곡선 (변별도 = 1 vs 2)
#   fig02a_crc_boundary        그림 2   범주반응곡선과 경계모수 위치
#   fig02b_boundary_curves     그림 2   인접 두 범주만 남긴 경계반응곡선
#   fig03_crc_transposition    그림 3   순서대로인 경계모수 vs 전치된 경계모수
#   fig04_scoring_function     그림 4   채점함수 간격 (수직선)
#   fig05_crc_scoring_function 그림 5   채점함수 간격에 따른 범주반응곡선
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[02] 연구 목적 개념 그림")

THETA <- seq(CONCEPT_XLIM[1], CONCEPT_XLIM[2], length.out = 401)

# 범주반응곡선 공통: 색 = 범주, 범주 번호를 곡선 꼭대기에 직접 표시 (범례 없음)
crc_layers <- function(df, group_vars = character(0)) {
  lab <- crc_peak_labels(df, group_vars)
  list(
    geom_line(aes(colour = category), linewidth = LW),
    geom_label(data = lab, aes(label = category, colour = category),
               vjust = -0.25, size = BASE_SIZE / 3.2, family = FONT_FAMILY,
               label.size = 0, fill = alpha("white", 0.8), show.legend = FALSE),
    scale_colour_manual(values = COL_CATEGORY, guide = "none"),
    scale_x_continuous(breaks = seq(CONCEPT_XLIM[1], CONCEPT_XLIM[2], 1)),
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)),
    labs(x = LAB$theta, y = LAB$prob),
    theme_thesis()
  )
}

b_line_layers <- function(b, y_text = 1.05) {
  bl <- data.frame(b = b, lab = paste0("b[", seq_along(b), "]"))
  list(
    geom_vline(data = bl, aes(xintercept = b), linetype = "dashed", colour = "grey35", linewidth = 0.8),
    geom_text(data = bl, aes(x = b, y = y_text, label = lab), parse = TRUE,
              size = BASE_SIZE / 3, family = FONT_FAMILY, hjust = -0.15)
  )
}

# ── 그림 1: 변별도에 따른 범주반응곡선 ───────────────────────────────────────
d1 <- bind_rows(lapply(CONCEPT_A_PAIR, function(a)
  crc_long(THETA, crc_true(THETA, CONCEPT_B, a = a), a = paste0("변별도 = ", a))))
d1$a <- factor(d1$a, levels = paste0("변별도 = ", CONCEPT_A_PAIR))
p1 <- ggplot(d1, aes(theta, prob)) +
  crc_layers(d1, "a") +
  facet_wrap(~a, nrow = 1)
save_fig(p1, "fig01_crc_discrimination", 14, 6.5)

# ── 그림 2a: 범주반응곡선 + 경계모수 ─────────────────────────────────────────
d2 <- crc_long(THETA, crc_true(THETA, CONCEPT_B))
p2a <- ggplot(d2, aes(theta, prob)) +
  b_line_layers(CONCEPT_B, y_text = 1.0) +
  crc_layers(d2)
save_fig(p2a, "fig02a_crc_boundary", 11, 7)

# ── 그림 2b: 경계반응곡선 P(범주 k+1 | 범주 k 또는 k+1) = logistic(a(θ − b_k)) ──
d2b <- bind_rows(lapply(1:4, function(k)
  data.frame(theta = THETA, prob = plogis(DISCRIM * (THETA - CONCEPT_B[k])),
             k = paste0("b", k))))
lab2b <- data.frame(k = paste0("b", 1:4), b = CONCEPT_B, lab = paste0("b[", 1:4, "]"))
p2b <- ggplot(d2b, aes(theta, prob, colour = k)) +
  geom_hline(yintercept = 0.5, colour = "grey50", linewidth = 0.6) +
  geom_vline(data = lab2b, aes(xintercept = b, colour = k), linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = LW) +
  geom_text(data = lab2b, aes(x = b, y = 1.05, label = lab), parse = TRUE,
            size = BASE_SIZE / 3, family = FONT_FAMILY, hjust = -0.15, show.legend = FALSE) +
  scale_colour_manual(values = COL_PARAM, guide = "none") +
  scale_x_continuous(breaks = seq(CONCEPT_XLIM[1], CONCEPT_XLIM[2], 1)) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(x = LAB$theta, y = "인접 두 범주 중 높은 범주의 응답 확률") +
  theme_thesis()
save_fig(p2b, "fig02b_boundary_curves", 11, 7)

# ── 그림 3: 순서대로인 경계모수 vs 전치된 경계모수 ───────────────────────────
panel3 <- c("b[3] < b[4]", "b[3] > b[4]")
d3 <- bind_rows(
  crc_long(THETA, crc_true(THETA, CONCEPT_B),            panel = panel3[1]),
  crc_long(THETA, crc_true(THETA, CONCEPT_B_TRANSPOSED), panel = panel3[2]))
d3$panel <- factor(d3$panel, levels = panel3)
bl3 <- bind_rows(
  data.frame(panel = panel3[1], b = CONCEPT_B,            lab = paste0("b[", 1:4, "]")),
  data.frame(panel = panel3[2], b = CONCEPT_B_TRANSPOSED, lab = paste0("b[", 1:4, "]")))
bl3$panel <- factor(bl3$panel, levels = panel3)
p3 <- ggplot(d3, aes(theta, prob)) +
  geom_vline(data = bl3, aes(xintercept = b), linetype = "dashed", colour = "grey35", linewidth = 0.8) +
  geom_text(data = bl3, aes(x = b, y = 1.0, label = lab), parse = TRUE,
            size = BASE_SIZE / 3.3, family = FONT_FAMILY, hjust = -0.1) +
  crc_layers(d3, "panel") +
  facet_wrap(~panel, nrow = 1, labeller = label_parsed)
save_fig(p3, "fig03_crc_transposition", 14, 6.5)

# ── 그림 4: 채점함수 간격 (수직선) ────────────────────────────────────────────
d4 <- bind_rows(lapply(sort(SIM1_IV1), function(s4)
  data.frame(sf4 = s4, k = 1:5, s = c(0, 1, 2, s4, 4))))
d4$row <- factor(unname(SF_GAP_LABEL[as.character(d4$sf4)]), levels = rev(SF_LEVELS))
seg4 <- d4 %>% filter(k %in% 4:5) %>% group_by(row) %>%
  summarise(x = min(s), xend = max(s), .groups = "drop")
p4 <- ggplot(d4, aes(s, row)) +
  geom_hline(aes(yintercept = as.numeric(row)), colour = "grey60", linewidth = 0.6) +
  geom_segment(data = seg4, aes(x = x, xend = xend, y = row, yend = row),
               colour = COL_CATEGORY[["5"]], linewidth = 4, alpha = 0.6) +
  geom_point(aes(colour = factor(k)), size = PS + 1.5) +
  geom_text(aes(label = paste0("s[", k, "]")), parse = TRUE, vjust = -1.1,
            size = BASE_SIZE / 3.2, family = FONT_FAMILY) +
  scale_colour_manual(values = COL_CATEGORY, guide = "none") +
  scale_x_continuous(breaks = 0:4, limits = c(-0.2, 4.2)) +
  scale_y_discrete(labels = lab_sf, expand = expansion(add = c(0.5, 0.9))) +
  labs(x = LAB$sf_value, y = NULL) +
  theme_thesis() +
  theme(panel.grid.major.y = element_blank())
save_fig(p4, "fig04_scoring_function", 12, 6)

# ── 그림 5: 채점함수 간격에 따른 범주반응곡선 ────────────────────────────────
d5 <- bind_rows(lapply(sort(SIM1_IV1), function(s4)
  crc_long(THETA, crc_true(THETA, CONCEPT_B, sf = c(0, 1, 2, s4, 4)), sf4 = s4)))
d5$f_sf <- factor(unname(SF_GAP_LABEL[as.character(d5$sf4)]), levels = SF_LEVELS)
p5 <- ggplot(d5, aes(theta, prob)) +
  crc_layers(d5, "f_sf") +
  facet_wrap(~f_sf, nrow = 1, labeller = as_labeller(lab_sf))
save_fig(p5, "fig05_crc_scoring_function", 16, 6.5)
