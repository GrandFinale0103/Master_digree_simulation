# =============================================================================
# 27_analysis_crossing_point.R  (simulation2, GPCM 추정 분석)
# 문항 1의 b3, b4 편향을 "교차점" 식으로 예측할 때,
# 교차점 대신 해당 범주 쌍 응답자의 평균 θ 를 쓰면 예측이 얼마나 달라지는지 비교
#
# ── 편향 예측식 ───────────────────────────────────────────────────────────────
#   생성 모형의 인접 범주 로짓 : w_k (θ − b_k)      w = 채점함수 간격 (b3: s4 − 2, b4: 4 − s4)
#   GPCM 극한의 인접 범주 로짓 : a (θ − b̂_k)       a, b̂ = 극한값 (17번 limit_values.csv)
#   두 직선이 만나는 θ = 교차점 θ*
#     θ* = (a · b̂_k − w_k · b_k) / (a − w_k)
#   편향 예측 = (w_k / a − 1) × (b_k − 기준점)
#     pred_star        : 기준점 = 교차점 θ*        (정의상 극한 편향 b̂_k − b_k 와 정확히 같음)
#     pred_mean        : 기준점 = 평균 θ (시뮬레이션) — 범주 k·k+1 응답자의 평균 θ (26번 집계)
#     pred_mean_theory : 기준점 = 평균 θ (이론값, θ ~ N(0,1) 적분)
#     pred_mean_a1     : 기준점 = 평균 θ (시뮬레이션), a = 1 로 둔 식 (w_k − 1)(b_k − 평균 θ)
#                        → 극한값 없이 응답 자료만으로 계산 가능한 예측
#   s4 = 3 조건은 w_k = a = 1 이라 교차점이 정의되지 않고 편향도 0 → 교차점 관련 값은 NA,
#   요약 지표에서 제외
#
# ── 비교 대상 (실제 편향) ─────────────────────────────────────────────────────
#   asym_bias : 극한 편향 b̂_k − b_k  (표본이 무한할 때)
#   obs_bias  : 시뮬레이션 GPCM 추정 편향 mean(추정값) − b_k  (17번 limit_bias_summary.csv 의 bias_true)
#
# ── 교차점과 평균 θ 의 차이 ───────────────────────────────────────────────────
#   diff_sim    = θ* − 평균 θ (시뮬레이션)
#   diff_theory = θ* − 평균 θ (이론값)
#
# 입력 (먼저 실행할 것):
#   output/analysis/limit/limit_values.csv          — 17_analysis_limit_bias.R (필수)
#   output/analysis/pair_side/pair_side_cond.csv    — 26_analysis_pair_side.R (필수, PAIRS 에 3·4 포함)
#   output/analysis/limit/limit_bias_summary.csv    — 17_analysis_limit_bias.R (없으면 obs_bias 비교 생략)
# 출력: output/analysis/crossing_point/
#   crossing_point_cond.csv          — 조건 × 경계(b3, b4)별 교차점, 평균 θ, 차이, 예측·실제 편향
#   crossing_point_summary.txt       — 예측 정확도 표, 차이 요약, 간격 × 난이도 표
#   crossing_point_diff_plot.png     — θ* − 평균 θ (행 IV1, 열 IV2 좁음→넓음, x축 IV3, 색 b3/b4)
#   crossing_point_values_plot_b3.png / _b4.png — 참값 b, 교차점, 평균 θ (같은 배치)
#   crossing_point_pred_scatter.png  — 예측 편향(x) 대 실제 편향(y), 기준점별
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

LIMIT_PATH <- "output/analysis/limit/limit_values.csv"
PAIR_PATH  <- "output/analysis/pair_side/pair_side_cond.csv"
OBS_PATH   <- "output/analysis/limit/limit_bias_summary.csv"
OUT_DIR    <- "output/analysis/crossing_point"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

if (!file.exists(LIMIT_PATH)) stop("먼저 R/17_analysis_limit_bias.R 을 실행하세요: ", LIMIT_PATH)
if (!file.exists(PAIR_PATH))  stop("먼저 R/26_analysis_pair_side.R 을 실행하세요: ", PAIR_PATH)

key_of <- function(d) sprintf("%.4f_%.4f_%.4f", d$iv1_sf4, d$iv2_b_interval, d$iv3_b_mean)

# =============================================================================
# 1. 입력 결합
# =============================================================================
limits <- read.csv(LIMIT_PATH)
limits$key <- key_of(limits)
lim_long <- limits %>%
  select(key, a_limit, b3_true, b4_true, b3_limit, b4_limit) %>%
  pivot_longer(-c(key, a_limit), names_to = c("param", ".value"), names_sep = "_")

pair <- read.csv(PAIR_PATH, colClasses = c(cond_code = "character")) %>%
  filter(boundary %in% 3:4) %>%
  mutate(key = key_of(.), param = paste0("b", boundary))
if (!all(c("b3", "b4") %in% pair$param))
  stop("26번 출력에 b3, b4 가 모두 있어야 합니다 (26번 PAIRS 에 3, 4 포함).")

df <- pair %>%
  select(cond_code, key, iv1_sf4, iv2_b_interval, iv3_b_mean, param, b_true,
         mean_theta, mean_theta_theory, prop_below_pooled) %>%
  left_join(lim_long, by = c("key", "param"))
if (any(is.na(df$a_limit))) stop("극한값과 맞지 않는 조건이 있습니다 (17번과 26번의 조건 확인).")
if (max(abs(df$b_true - df$true)) > 1e-6) stop("17번과 26번의 참 경계모수가 다릅니다.")

if (file.exists(OBS_PATH)) {
  obs <- read.csv(OBS_PATH) %>% mutate(key = key_of(.)) %>%
    select(key, param, obs_bias = bias_true, n_obs = n)
  df <- df %>% left_join(obs, by = c("key", "param"))
} else {
  message("참고: ", OBS_PATH, " 없음 → 시뮬레이션 편향(obs_bias) 비교 생략")
  df$obs_bias <- NA_real_; df$n_obs <- NA_integer_
}

# =============================================================================
# 2. 교차점, 차이, 예측 편향
# =============================================================================
df <- df %>%
  mutate(w         = ifelse(param == "b3", iv1_sf4 - 2, 4 - iv1_sf4),
         a         = a_limit,
         b_limit   = limit,
         asym_bias = b_limit - b_true,
         defined   = abs(a - w) > 1e-3,                      # s4 = 3 이면 FALSE
         theta_star       = ifelse(defined, (a * b_limit - w * b_true) / (a - w), NA_real_),
         diff_sim         = theta_star - mean_theta,
         diff_theory      = theta_star - mean_theta_theory,
         pred_star        = ifelse(defined, (w / a - 1) * (b_true - theta_star), NA_real_),
         pred_mean        = (w / a - 1) * (b_true - mean_theta),
         pred_mean_theory = (w / a - 1) * (b_true - mean_theta_theory),
         pred_mean_a1     = (w - 1)     * (b_true - mean_theta)) %>%
  select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param, b_true, w, a, b_limit,
         theta_star, mean_theta, mean_theta_theory, diff_sim, diff_theory, prop_below_pooled,
         pred_star, pred_mean, pred_mean_theory, pred_mean_a1, asym_bias, obs_bias, n_obs, defined) %>%
  arrange(cond_code, param)
write.csv(df, out("crossing_point_cond.csv"), row.names = FALSE)

# =============================================================================
# 3. 요약
# =============================================================================
PREDS <- c(pred_star        = "교차점 θ*",
           pred_mean        = "평균 θ (시뮬레이션)",
           pred_mean_theory = "평균 θ (이론값)",
           pred_mean_a1     = "평균 θ (시뮬레이션), a = 1")

acc <- function(pred, act) {
  ok <- is.finite(pred) & is.finite(act)
  p <- pred[ok]; y <- act[ok]
  if (length(y) < 3) return(c(n = length(y), cor = NA, R2 = NA, slope = NA, MAE = NA, sign_agree = NA))
  c(n = length(y), cor = cor(p, y), R2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2),
    slope = unname(coef(lm(y ~ 0 + p))), MAE = mean(abs(y - p)), sign_agree = mean(sign(p) == sign(y)))
}

acc_tab <- df %>% filter(defined) %>%
  pivot_longer(names(PREDS), names_to = "predictor", values_to = "pred") %>%
  pivot_longer(c(asym_bias, obs_bias), names_to = "target", values_to = "actual") %>%
  group_by(param, predictor, target) %>%
  summarise(m = list(acc(pred, actual)), .groups = "drop") %>%
  unnest_wider(m) %>%
  mutate(predictor = factor(predictor, levels = names(PREDS)),
         target    = factor(target, levels = c("asym_bias", "obs_bias"))) %>%
  arrange(param, target, predictor)
write.csv(acc_tab, out("crossing_point_accuracy.csv"), row.names = FALSE)

diff_tab <- df %>% filter(defined) %>%
  group_by(param) %>%
  summarise(n = n(),
            mean_diff_sim = mean(diff_sim), mean_abs_diff_sim = mean(abs(diff_sim)),
            min_diff_sim = min(diff_sim), max_diff_sim = max(diff_sim),
            mean_diff_theory = mean(diff_theory), mean_abs_diff_theory = mean(abs(diff_theory)),
            cor_star_mean = cor(theta_star, mean_theta), .groups = "drop")

sink(out("crossing_point_summary.txt"), split = TRUE)
cat("교차점 θ* 대 평균 θ — simulation2 문항 1 (s4 = 3 조건 제외)\n")
cat("편향 예측 = (w/a − 1)(b − 기준점), a = 1 식만 (w − 1)(b − 평균 θ)\n")
cat("pred_star 는 정의상 극한 편향(asym_bias)과 같음 — 기준선으로 표시\n")
cat("\n[1] θ* − 평균 θ (양수 = 교차점이 평균 θ 보다 위)\n")
print(as.data.frame(diff_tab %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

for (p in c("b3", "b4")) {
  cat(sprintf("\n  %s: θ* − 평균 θ (시뮬레이션), 간격(행) × 난이도(열), s4 수준 평균\n", p))
  d <- df %>% filter(defined, param == p)
  print(round(tapply(d$diff_sim, list(interval = d$iv2_b_interval, mean = d$iv3_b_mean), mean), 2))
}

cat("\n[2] 편향 예측 정확도\n")
cat("  target asym_bias = 극한 편향, obs_bias = 시뮬레이션 GPCM 편향\n")
cat("  R2 = 1 − SSE/SST (예측값을 그대로 쓸 때), slope = 실제 ~ 0 + 예측 기울기 (1이면 크기 일치)\n")
cat("  sign_agree = 예측과 실제의 부호가 같은 조건 비율\n")
cat("  predictor: pred_star = 교차점, pred_mean = 평균 θ (시뮬레이션), pred_mean_theory = 평균 θ (이론값),\n")
cat("             pred_mean_a1 = 평균 θ (시뮬레이션) + a = 1\n")
old_width <- options(width = 200)
print(as.data.frame(acc_tab %>%
        mutate(across(c(cor, R2, slope, MAE, sign_agree), ~ round(.x, 3)))), row.names = FALSE)
options(old_width)
sink()

# =============================================================================
# 4. 그래프
# =============================================================================
make_num_factor <- function(x, fmt) {
  vals <- sort(unique(na.omit(x)))
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
grid_layers <- list(
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                 iv2 = function(x) paste0("Int=", x))),
  labs(x = "b-parameter Mean (IV3)"),
  theme_bw(base_size = 9),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
plot_df <- df %>% filter(defined) %>%
  mutate(iv1 = make_num_factor(iv1_sf4, "%.2f"),
         iv2 = make_num_factor(iv2_b_interval, "%.2f"))

# 4-1. θ* − 평균 θ
p_diff <- ggplot(plot_df, aes(iv3_b_mean, diff_sim, colour = param)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  geom_line(linewidth = 0.4) + geom_point(size = 1.2) +
  scale_colour_manual(values = c(b3 = "#F2A541", b4 = "#D1495B"), name = "Boundary") +
  scale_x_continuous(breaks = sort(unique(df$iv3_b_mean))) +
  labs(y = "crossing point - mean theta of respondents in the category pair",
       title = "Crossing point (theta*) minus mean theta (simulation), item 1") +
  grid_layers
ggsave(out("crossing_point_diff_plot.png"), p_diff, width = 24, height = 18, dpi = 150)

# 4-2. 참값 b, 교차점, 평균 θ
for (p in c("b3", "b4")) {
  v <- plot_df %>% filter(param == p) %>%
    select(iv1, iv2, iv3_b_mean, b_true, theta_star, mean_theta) %>%
    pivot_longer(c(b_true, theta_star, mean_theta), names_to = "what", values_to = "value") %>%
    mutate(what = factor(what, levels = c("b_true", "theta_star", "mean_theta"),
                         labels = c(sprintf("true %s", p), "crossing point theta*", "mean theta (pair)")))
  pv <- ggplot(v, aes(iv3_b_mean, value, colour = what, linetype = what)) +
    geom_line(linewidth = 0.4) + geom_point(size = 1.1) +
    scale_colour_manual(values = c("grey30", "#2a6fb0", "#c0392b"), name = NULL) +
    scale_linetype_manual(values = c("dashed", "solid", "solid"), name = NULL) +
    scale_x_continuous(breaks = sort(unique(df$iv3_b_mean))) +
    labs(y = "theta", title = sprintf("Item 1 %s: true value, crossing point and mean theta", p)) +
    grid_layers
  ggsave(out(sprintf("crossing_point_values_plot_%s.png", p)), pv, width = 24, height = 18, dpi = 150)
}

# 4-3. 예측 편향 대 실제 편향
target_col <- if (all(is.na(df$obs_bias))) "asym_bias" else "obs_bias"
target_lab <- if (target_col == "obs_bias") "simulation GPCM bias (mean estimate - true)" else
              "asymptotic bias (limit - true)"
sc <- plot_df %>%
  pivot_longer(names(PREDS), names_to = "predictor", values_to = "pred") %>%
  mutate(actual = .data[[target_col]],
         predictor = factor(predictor, levels = names(PREDS),
                            labels = c("crossing point", "mean theta (sim)",
                                       "mean theta (theory)", "mean theta (sim), a = 1")))
p_sc <- ggplot(sc, aes(pred, actual, colour = factor(iv1_sf4))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 1.3, alpha = 0.8) +
  facet_grid(param ~ predictor) +
  scale_colour_viridis_d(name = "s4", end = 0.9) +
  labs(x = "predicted bias = (w/a - 1)(b - reference point)", y = target_lab,
       title = "Predicted vs actual boundary bias by reference point (dashed = y = x)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(out("crossing_point_pred_scatter.png"), p_sc, width = 14, height = 8, dpi = 150)

message(sprintf("완료: %s", OUT_DIR))
