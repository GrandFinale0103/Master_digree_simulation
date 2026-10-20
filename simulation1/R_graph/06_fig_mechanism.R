# =============================================================================
# 06_fig_mechanism.R  (simulation1/R_graph)
# 결과 — GPCM 편향이 생기는 이유 (시뮬레이션1)
#
#   fig10_estimate_distribution   그림 10    b̂3·b̂4 추정값 분포와 참값 (FIG10_CONDS 조건들, 같은 x축)
#   fig11_1_b4_curves             그림 11-1  범주 4·5 응답자의 참 / GPCM 문항특성곡선 (b4)
#   fig11_2_b4_curves_theta       그림 11-2  11-1 + 범주 4·5 응답자 θ (초록: θ < b4, 파랑: θ ≥ b4)
#   fig11_3_b3_curves             그림 11-3  범주 3·4 응답자의 참 / GPCM 문항특성곡선 (b3)
#   fig11_4_b3_curves_theta       그림 11-4  11-3 + 범주 3·4 응답자 θ (초록: θ < b3, 파랑: θ ≥ b3)
#   extra_kl_vs_sim_bias          예비       KL 극한으로 계산한 Bias vs 시뮬레이션 Bias
#
# 그림 11 의 곡선 (두 범주만 남겼을 때 높은 범주를 고를 확률)
#   참 모형  : logistic(a · (s_{k+1} − s_k) · (θ − b_k)),  a = 1, 채점함수 (0, 1, 2, s4, 4)
#   GPCM     : logistic(â · (θ − b_k)),  â = 그 조건 반복들의 GPCM 변별도 평균
#              위치는 b̂ 대신 참 b 에 맞춤 → 두 곡선의 기울기 차이만 남아
#              b 아래·위에서 GPCM 이 logit 을 과대·과소 예측하는 것을 설명할 수 있음
#              (실제 추정 평균 b̂ 는 tab_fig11_summary 에 기록)
#   점       : 대표 반복 1개의 응답 파일에서 문항 1에 두 범주 중 하나로 응답한 사람의 참 θ
#              (FIG11_REP = NULL 이면 b̂3·b̂4 가 조건 평균에 가장 가까운 반복, 유령 응답자 제외)
#              모두 y = 0 에 찍고 투명도(FIG11_POINT_ALPHA)를 줘서 겹칠수록 진해짐
#   네 그림 모두 x축·y축 범위가 같음
#
# 표
#   tab_fig11_summary   그림 11 에 쓴 값 (참값, 추정 평균, 대표 반복, 초록·파랑 점 수)
#   tab_kl_limit        조건 × 경계모수별 참값, KL 극한, 시뮬레이션 평균
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")
message("[06] 편향 원리")

reps1 <- prepare_reps(load_sim1())

# =============================================================================
# 그림 10: 추정값 분포
# =============================================================================
d10 <- bind_rows(lapply(FIG10_CONDS, function(cd) {
  reps1[match_cond(reps1, cd), ] %>%
    transmute(panel = cond_text(cd), b3, b4, b3_true, b4_true)
}))
if (nrow(d10) == 0) {
  message("  그림 10: FIG10_CONDS 조건의 자료가 없음 → 건너뜀")
} else {
  d10$panel <- factor(d10$panel, levels = unique(d10$panel))
  est10  <- d10 %>% pivot_longer(c(b3, b4), names_to = "param", values_to = "est")
  true10 <- d10 %>% distinct(panel, b3 = b3_true, b4 = b4_true) %>%
    pivot_longer(c(b3, b4), names_to = "param", values_to = "true")
  hat_lab <- c(b3 = expression(hat(b)[3]), b4 = expression(hat(b)[4]))
  p10 <- ggplot(est10, aes(est, fill = param, colour = param)) +
    geom_density(alpha = 0.35, linewidth = LW, adjust = 1.2) +
    geom_vline(data = true10, aes(xintercept = true, colour = param, linetype = "true"),
               linewidth = LW, key_glyph = "path",
               show.legend = c(colour = FALSE, fill = FALSE, linetype = TRUE)) +
    scale_fill_manual(values = COL_PARAM, labels = hat_lab, name = NULL) +
    scale_colour_manual(values = COL_PARAM, labels = hat_lab, name = NULL) +
    scale_linetype_manual(values = c(true = "dashed"), labels = "참값", name = NULL) +
    facet_wrap(~panel, ncol = 1) +
    labs(x = LAB$estimate, y = LAB$density) +
    theme_thesis() +
    guides(linetype = guide_legend(override.aes = list(colour = "black")))
  save_fig(p10, "fig10_estimate_distribution", 15, 4 + 4 * nlevels(d10$panel))
}

# =============================================================================
# 그림 11: 문항특성곡선과 응답자 θ
# =============================================================================
pick_rep <- function(cr) {
  if (!is.null(FIG11_REP)) return(FIG11_REP)
  z <- ((cr$b3 - mean(cr$b3)) / sd(cr$b3))^2 + ((cr$b4 - mean(cr$b4)) / sd(cr$b4))^2
  cr$rep_id[which.min(z)]
}

read_response <- function(cond, rep) {
  code  <- sim1_code(cond[["sf4"]], cond[["interval"]], cond[["mean"]])
  files <- list.files(file.path(SIM1_OUTPUT, "responses"),
                      pattern = sprintf("_cond%s_rep%04d_.*_response\\.csv$", code, rep),
                      full.names = TRUE)
  if (length(files) == 0) stop(sprintf("응답 파일이 없습니다: 조건 %s, 반복 %d (%s)",
                                       code, rep, file.path(SIM1_OUTPUT, "responses")))
  read.csv(sort(files)[length(files)])
}

THETA11 <- seq(FIG11_XLIM[1], FIG11_XLIM[2], length.out = 400)
Y11     <- c(0, 1.06)
fig11_rows <- list()

make_fig11 <- function(cond, k, fig_a, fig_b) {
  cr <- reps1[match_cond(reps1, cond), ]
  if (nrow(cr) == 0) { message("  그림 11: 조건 자료 없음 → 건너뜀 (", cond_text(cond), ")"); return() }
  sf     <- c(0, 1, 2, cond[["sf4"]], 4)
  w_k    <- sf[k + 1] - sf[k]
  b_true <- cr[[paste0("b", k, "_true")]][1]
  a_hat  <- mean(cr$a, na.rm = TRUE)
  b_hat  <- mean(cr[[paste0("b", k)]])

  # 두 곡선 모두 참 b 에서 0.5 를 지나게 그림 → 기울기 차이만 남아
  # b 보다 아래·위에서 GPCM 이 logit 을 과대·과소 예측하는 모습을 보여 줌
  curves <- bind_rows(
    data.frame(theta = THETA11, prob = plogis(DISCRIM * w_k * (THETA11 - b_true)), model = "true"),
    data.frame(theta = THETA11, prob = plogis(a_hat * (THETA11 - b_true)),         model = "gpcm"))
  curves$model <- factor(curves$model, levels = c("true", "gpcm"))

  base <- ggplot() +
    geom_hline(yintercept = c(0, 0.5), colour = "grey60", linewidth = 0.5) +
    geom_vline(xintercept = b_true, colour = "grey30", linetype = "dashed", linewidth = 0.9) +
    annotate("text", x = b_true, y = 1.04, label = sprintf("b[%d] == hat(b)[%d]", k, k), parse = TRUE,
             hjust = -0.08, size = BASE_SIZE / 3, family = FONT_FAMILY) +
    geom_line(data = curves, aes(theta, prob, colour = model, linetype = model), linewidth = LW + 0.3) +
    scale_colour_manual(values = c(true = COL_TRUE, gpcm = COL_GPCM),
                        labels = c(true = LAB$true_model, gpcm = LAB$gpcm_fig11), name = NULL) +
    scale_linetype_manual(values = c(true = "solid", gpcm = "longdash"),
                          labels = c(true = LAB$true_model, gpcm = LAB$gpcm_fig11), name = NULL) +
    scale_x_continuous(breaks = seq(ceiling(FIG11_XLIM[1]), floor(FIG11_XLIM[2]), 1),
                       limits = FIG11_XLIM) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = Y11) +
    labs(x = LAB$theta, y = LAB$prob) +
    theme_thesis()
  save_fig(base, fig_a, 12, 8.5)   # 점 그림과 같은 크기 → 슬라이드에서 겹쳐도 곡선 위치가 같음

  # 점: 대표 반복에서 범주 k 또는 k+1 로 응답한 사람 (응답값 0~4 → 범주 = 응답값 + 1)
  #     모두 y = 0 에 찍고 투명도를 줘서 겹칠수록 진해짐
  rep  <- pick_rep(cr)
  resp <- read_response(cond, rep)
  pts  <- resp %>% filter(!is.na(theta), item1 %in% c(k - 1, k)) %>%
    mutate(side = factor(ifelse(theta < b_true, "below", "above"), levels = c("below", "above")))
  n_out <- sum(pts$theta < FIG11_XLIM[1] | pts$theta > FIG11_XLIM[2])
  if (n_out > 0) message(sprintf("  [주의] %s: x축 범위 밖 점 %d개", fig_b, n_out))

  side_lab <- c(below = parse(text = sprintf("theta < b[%d]", k))[[1]],
                above = parse(text = sprintf("theta >= b[%d]", k))[[1]])
  p_b <- base +
    geom_point(data = pts, aes(theta, 0, fill = side), shape = 21, colour = "transparent",
               size = PS + 0.5, alpha = FIG11_POINT_ALPHA) +
    scale_fill_manual(values = c(below = COL_BELOW, above = COL_ABOVE),
                      labels = side_lab, name = NULL, drop = FALSE) +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 1),
           fill = guide_legend(order = 2, override.aes = list(size = PS + 1, alpha = 1)))
  save_fig(p_b, fig_b, 12, 8.5)

  fig11_rows[[length(fig11_rows) + 1]] <<- data.frame(
    `그림` = paste(fig_a, fig_b, sep = ", "), `경계모수` = paste0("b", k),
    `조건` = cond_text(cond, sep = ", "),
    `참값` = round(b_true, 3), `참 기울기` = round(DISCRIM * w_k, 3),
    `GPCM 추정 평균` = round(b_hat, 3), `GPCM 변별도 평균` = round(a_hat, 3),
    `대표 반복` = rep, `응답자 수` = nrow(pts),
    `초록(θ < b)` = sum(pts$side == "below"), `파랑(θ ≥ b)` = sum(pts$side == "above"),
    check.names = FALSE)
}

make_fig11(FIG11_B4_COND, 4, "fig11_1_b4_curves", "fig11_2_b4_curves_theta")
make_fig11(FIG11_B3_COND, 3, "fig11_3_b3_curves", "fig11_4_b3_curves_theta")
if (length(fig11_rows)) save_table(bind_rows(fig11_rows), "tab_fig11_summary")

# =============================================================================
# 예비: KL 극한 편향 vs 시뮬레이션 편향
#   모형이 틀렸을 때 최대우도 추정값이 수렴하는 값 (θ ~ N(0, 1), simulation1/R/13번과 같은 계산)
# =============================================================================
Q  <- seq(-6, 6, length.out = 1201)
WQ <- dnorm(Q) / sum(dnorm(Q))
kl_limit <- function(sf4, b) {
  P   <- crc_true(Q, b, sf = c(0, 1, 2, sf4, 4))
  obj <- function(par) -sum(WQ * rowSums(P * log(pmax(crc_true(Q, par[2:5], a = exp(par[1])), 1e-300))))
  o   <- optim(c(0, b), obj, method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
  setNames(o$par[2:5], paste0("b", 1:4))
}

conds <- reps1 %>% distinct(sf4, interval, mean, f_sf, f_int, f_mean)
kl <- bind_rows(lapply(seq_len(nrow(conds)), function(i) {
  cd  <- conds[i, ]
  b   <- cd$mean + cd$interval * B_OFFSET
  lim <- kl_limit(cd$sf4, b)
  sim <- reps1[match_cond(reps1, c(sf4 = cd$sf4, interval = cd$interval, mean = cd$mean)), ]
  data.frame(cd, param = paste0("b", 1:4), true = b, limit = unname(lim),
             sim_mean = sapply(paste0("b", 1:4), function(p) mean(sim[[p]])), row.names = NULL)
})) %>% mutate(kl_bias = limit - true, sim_bias = sim_mean - true)

rng <- range(c(kl$kl_bias, kl$sim_bias))
param_lab <- c(b1 = expression(b[1]), b2 = expression(b[2]), b3 = expression(b[3]), b4 = expression(b[4]))
p_kl <- ggplot(kl, aes(kl_bias, sim_bias, colour = param, shape = param)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linewidth = 0.8) +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = PS + 0.5, alpha = 0.9) +
  scale_colour_manual(values = COL_PARAM, labels = param_lab, name = NULL) +
  scale_shape_manual(values = c(b1 = 16, b2 = 17, b3 = 15, b4 = 18), labels = param_lab, name = NULL) +
  coord_equal(xlim = rng, ylim = rng) +
  labs(x = LAB$kl_bias, y = LAB$sim_bias) +
  theme_thesis()
save_fig(p_kl, "extra_kl_vs_sim_bias", 10, 10)

save_table(kl %>% transmute(`채점함수 간격` = f_sf, `경계모수 간격` = f_int, `문항 난이도` = f_mean,
                            `경계모수` = param, `참값` = round(true, 3), `KL 극한` = round(limit, 3),
                            `추정 평균` = round(sim_mean, 3), `KL 극한 Bias` = round(kl_bias, 3),
                            `시뮬레이션 Bias` = round(sim_bias, 3)),
           "tab_kl_limit")
