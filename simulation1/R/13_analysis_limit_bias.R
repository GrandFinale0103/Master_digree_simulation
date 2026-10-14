# =============================================================================
# 13_analysis_limit_bias.R  (simulation1, GPCM 추정 결과)
# 채점함수 오지정 GPCM 추정에서 문항 1의 b3, b4 가 수렴하는 값(극한값)을 계산하고,
# 그 값을 기준으로 Bias / RMSE 를 계산·시각화
#
# ── 극한값(pseudo-true value) ────────────────────────────────────────────────
#   모형이 틀려도 최대우도 추정값은 참값이 아니라
#   "참 모형과 추정 모형의 KL 발산을 최소로 만드는 모수"로 수렴함 (White, 1982)
#     참 모형  : 채점함수 (0, 1, 2, s4, 4), a = 1, b = 평균 + 간격 × (−1.5, −0.5, 0.5, 1.5)
#     추정 모형: 채점함수 (0, 1, 2, 3, 4) GPCM (a, b1~b4)
#   극한값 = argmax_{a,b} ∫ g(θ) Σ_k P_true(k|θ) log P_GPCM(k|θ; a, b) dθ
#     g(θ) = 조건의 능력모수 분포 (정규, 정적·부적편포, 균등)
#   (문항 2~20 은 올바르게 지정되어 θ 척도를 정하므로 문항 1만 따로 계산)
#
# ── 척도 (균등분포 조건 주의) ─────────────────────────────────────────────────
#   mirt 는 θ ~ N(0, 1) 로 척도를 고정하므로 추정값은 평균 0, 분산 1 척도에 있음
#   정규·편포 조건은 생성 분포도 평균 0, 분산 1 → 생성 척도 = 추정 척도
#   균등분포 U(−4, 4) 는 표준편차 σ = 8/√12 ≈ 2.31 → 추정 척도의 참값은 b/σ, a·σ
#     b_true     : 생성 척도 참값 (09번 Bias 의 기준)
#     b_true_std : 추정 척도 참값 (균등분포만 b/σ, 나머지는 b_true 와 같음)
#     극한값은 추정 척도에서 계산 (θ 분포를 평균 0, 분산 1 로 표준화)
#   정규가 아닌 분포는 mirt 의 정규 사전분포와 모양이 달라 "kl" 극한값이 근사임
#   → 정확한 값이 필요하면 LIMIT_METHOD = "mirt_large" 로 확인
#
#   LIMIT_METHOD
#     "kl"         : 위 적분을 구적법으로 계산 (조건당 1초 안팎, 기본값)
#     "mirt_large" : 조건마다 N_LARGE 명의 응답을 생성해 mirt GPCM 으로 한 번 추정한 값
#                    (실제 추정 절차와 같음, 조건당 수십 초 — 결과를 저장해 두고 재사용)
#
# ── 지표 (수렴한 반복, 조건별) ────────────────────────────────────────────────
#   bias_limit = mean(b̂ − b_limit)            — 극한값 기준 Bias (유한 표본 편향)
#   rmse_limit = sqrt(mean((b̂ − b_limit)²))   — 극한값 기준 RMSE
#   asym_bias  = b_limit − b_true_std          — 표본이 무한히 커도 남는 Bias (채점함수 오지정 몫)
#   bias_true / rmse_true                      — 생성 척도 참값 기준 (09번과 같은 값, 비교용)
#
# 입력: output/analysis/rmse_bias_long.csv  (09_analysis_rmse_bias.R 출력)
# 출력: output/analysis/limit/
#   limit_values.csv          — 108 조건 전체의 극한값 (a, b1~b4) 과 참값
#   limit_bias_summary.csv    — 조건 × 모수(b3, b4)별 지표
#   limit_bias_plot.png       — 극한값 기준 Bias  (행 모수, 열 IV1 × IV2, x축 IV3, 색 IV4)
#   limit_rmse_plot.png       — 극한값 기준 RMSE  (같은 배치)
#   limit_asym_bias_plot.png  — 극한 Bias (b_limit − b_true_std)  (같은 배치)
#   limit_values_plot_b3.png / _b4.png — 참값·극한값·추정 평균 비교 (행 IV4, 열 IV1 × IV2)
# 정규분포 조건만: output/analysis_normal/limit/  (11_analysis_normal_only.R 과 같은 배치, 4단계 참고)
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
LIMIT_METHOD <- "kl"     # "kl" 또는 "mirt_large"
N_LARGE      <- 100000   # mirt_large: 조건당 생성할 응답자 수
LARGE_SEED   <- 20261014 # mirt_large: 난수 시드
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

IN_PATH <- "output/analysis/rmse_bias_long.csv"
OUT_DIR <- "output/analysis/limit"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

IV1_LEVELS <- c(3, 3.33, 3.66)
IV2_LEVELS <- c(1.5, 1, 0.5)
IV3_LEVELS <- c(0, 0.5, 1)
IV4_LEVELS <- c("pos_skew", "normal", "neg_skew", "uniform")
IV4_LABEL  <- c(pos_skew = "정적편포(+0.8)", normal = "정규분포",
                neg_skew = "부적편포(-0.8)", uniform = "균등분포")
B_OFFSET   <- c(-1.5, -0.5, 0.5, 1.5)
UNIF_SD    <- 8 / sqrt(12)                      # U(−4, 4) 의 표준편차

# ── 편포 분포 모수 (00_setup.R 의 find_rsn_params 와 같은 계산) ──────────────
find_rsn_params <- function(target_skew) {
  f_alpha <- function(alpha) {
    delta <- alpha / sqrt(1 + alpha^2)
    ((4 - pi) / 2) * (delta * sqrt(2 / pi))^3 / (1 - 2 * delta^2 / pi)^(3 / 2) - target_skew
  }
  interval <- if (target_skew > 0) c(0.001, 50) else c(-50, -0.001)
  alpha <- uniroot(f_alpha, interval = interval, tol = 1e-10)$root
  delta <- alpha / sqrt(1 + alpha^2)
  omega <- 1 / sqrt(1 - 2 * delta^2 / pi)
  list(alpha = alpha, xi = -omega * delta * sqrt(2 / pi), omega = omega)
}
RSN_POS <- find_rsn_params(0.8)
RSN_NEG <- find_rsn_params(-0.8)
dsn_std <- function(x, p) {                     # 왜도 정규분포 밀도
  z <- (x - p$xi) / p$omega
  2 / p$omega * dnorm(z) * pnorm(p$alpha * z)
}

# 추정 척도(평균 0, 분산 1)의 θ 밀도
theta_density <- function(x, dist) {
  switch(dist,
    normal   = dnorm(x),
    pos_skew = dsn_std(x, RSN_POS),
    neg_skew = dsn_std(x, RSN_NEG),
    uniform  = as.numeric(abs(x) <= sqrt(3)) / (2 * sqrt(3)),
    stop("알 수 없는 분포: ", dist))
}
# 추정 척도로 옮긴 문항 1 참값 (균등분포만 척도 변환)
scale_sd <- function(dist) if (dist == "uniform") UNIF_SD else 1

# =============================================================================
# 1. 극한값 계산
# =============================================================================
# eta_k = slope_k · θ − int_k 인 범주 확률 (행: θ, 열: 범주 0~4)
cat_probs <- function(theta, slope, int) {
  z <- outer(theta, slope) - matrix(int, length(theta), length(slope), byrow = TRUE)
  z <- z - apply(z, 1, max)
  p <- exp(z)
  p / rowSums(p)
}
# 참 모형: 채점함수 sf 의 인접 차이 w 를 쓰는 GPCM (03_response.R 과 같은 식)
true_probs <- function(theta, sf, b, a = 1) {
  w <- diff(sf)
  cat_probs(theta, a * c(0, cumsum(w)), a * c(0, cumsum(w * b)))
}
# 추정 모형: 표준 채점함수 GPCM, par = (log a, b1, b2, b3, b4)
gpcm_probs <- function(theta, par) {
  a <- exp(par[1])
  cat_probs(theta, a * 0:4, a * c(0, cumsum(par[2:5])))
}

QUAD <- seq(-6, 6, length.out = 1201)

limit_kl <- function(sf4, b, dist) {
  sd_g <- scale_sd(dist)
  w    <- theta_density(QUAD, dist); w <- w / sum(w)
  keep <- w > 0
  q    <- QUAD[keep]; w <- w[keep]
  b_s  <- b / sd_g                                 # 추정 척도 참값
  P    <- true_probs(q, c(0, 1, 2, sf4, 4), b_s, a = sd_g)
  obj  <- function(par) -sum(w * rowSums(P * log(pmax(gpcm_probs(q, par), 1e-300))))
  o    <- optim(c(log(sd_g), b_s), obj, method = "BFGS",
                control = list(reltol = 1e-14, maxit = 2000))
  if (o$convergence != 0)
    warning(sprintf("KL 최적화 미수렴: sf4=%g, b=%s, %s", sf4, paste(b, collapse = ","), dist))
  c(a = exp(o$par[1]), b1 = o$par[2], b2 = o$par[3], b3 = o$par[4], b4 = o$par[5])
}

limit_mirt_large <- function(sf4, b, dist, seed) {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("mirt 패키지가 필요합니다.")
  if (dist %in% c("pos_skew", "neg_skew") && !requireNamespace("sn", quietly = TRUE))
    stop("sn 패키지가 필요합니다.")
  set.seed(seed)
  theta <- switch(dist,                            # 03_response.R 의 sample_theta 와 같음
    normal   = rnorm(N_LARGE),
    pos_skew = sn::rsn(N_LARGE, xi = RSN_POS$xi, omega = RSN_POS$omega, alpha = RSN_POS$alpha),
    neg_skew = sn::rsn(N_LARGE, xi = RSN_NEG$xi, omega = RSN_NEG$omega, alpha = RSN_NEG$alpha),
    uniform  = runif(N_LARGE, -4, 4))
  theta  <- as.numeric(theta)
  item_b <- c(list(b), lapply(rnorm(19), function(m) m + 1.5 * B_OFFSET))
  resp <- sapply(seq_along(item_b), function(j) {
    sf <- if (j == 1) c(0, 1, 2, sf4, 4) else 0:4
    P  <- true_probs(theta, sf, item_b[[j]])
    1L + rowSums(runif(N_LARGE) > t(apply(P, 1, cumsum))[, 1:4])   # 범주 1~5
  })
  colnames(resp) <- paste0("item", 1:20)
  mod <- mirt::mirt(as.data.frame(resp), 1, itemtype = "gpcm", verbose = FALSE,
                    technical = list(NCYCLES = 3000))
  co <- mirt::coef(mod, IRTpars = TRUE, simplify = TRUE)$items[1, ]
  c(a = unname(co["a"]), b1 = unname(co["b1"]), b2 = unname(co["b2"]),
    b3 = unname(co["b3"]), b4 = unname(co["b4"]))
}

cond_grid <- expand.grid(iv4_theta_dist = IV4_LEVELS, iv3_b_mean = IV3_LEVELS,
                         iv2_b_interval = IV2_LEVELS, iv1_sf4 = IV1_LEVELS,
                         stringsAsFactors = FALSE)[, 4:1]
key <- function(d) sprintf("%.4f_%.4f_%.4f_%s", d$iv1_sf4, d$iv2_b_interval,
                           d$iv3_b_mean, d$iv4_theta_dist)
limit_path <- out(sprintf("limit_values_%s.csv", LIMIT_METHOD))
done <- if (file.exists(limit_path)) read.csv(limit_path, stringsAsFactors = FALSE) else NULL
todo <- cond_grid[!key(cond_grid) %in% if (is.null(done)) character(0) else key(done), ]

message(sprintf("극한값 계산 (%s): %d 조건 (저장된 값 %d 조건 재사용)",
                LIMIT_METHOD, nrow(todo), nrow(cond_grid) - nrow(todo)))
for (i in seq_len(nrow(todo))) {
  cg <- todo[i, ]
  b  <- cg$iv3_b_mean + cg$iv2_b_interval * B_OFFSET
  lim <- if (LIMIT_METHOD == "kl") limit_kl(cg$iv1_sf4, b, cg$iv4_theta_dist) else
         limit_mirt_large(cg$iv1_sf4, b, cg$iv4_theta_dist, LARGE_SEED + i)
  sd_g <- scale_sd(cg$iv4_theta_dist)
  row <- data.frame(cg, a_true = 1,
                    b1_true = b[1], b2_true = b[2], b3_true = b[3], b4_true = b[4],
                    a_true_std = sd_g, b3_true_std = b[3] / sd_g, b4_true_std = b[4] / sd_g,
                    a_limit = lim[["a"]], b1_limit = lim[["b1"]], b2_limit = lim[["b2"]],
                    b3_limit = lim[["b3"]], b4_limit = lim[["b4"]],
                    stringsAsFactors = FALSE)
  done <- bind_rows(done, row)
  if (LIMIT_METHOD == "mirt_large" || i == nrow(todo))
    write.csv(done, limit_path, row.names = FALSE)   # mirt_large 는 조건마다 저장 (중단 후 이어서)
  if (LIMIT_METHOD == "mirt_large")
    message(sprintf("  %d/%d  sf4=%g 간격=%g 평균=%g %s", i, nrow(todo),
                    cg$iv1_sf4, cg$iv2_b_interval, cg$iv3_b_mean, cg$iv4_theta_dist))
}
limits <- done %>% mutate(gap34_limit = b4_limit - b3_limit)
write.csv(limits, out("limit_values.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s (%d 조건, 극한 b4 − b3 < 0 인 조건 %d개)",
                out("limit_values.csv"), nrow(limits), sum(limits$gap34_limit < 0)))

# =============================================================================
# 2. 극한값 기준 Bias / RMSE
# =============================================================================
if (!file.exists(IN_PATH)) stop("먼저 R/09_analysis_rmse_bias.R 을 실행하세요: ", IN_PATH)
est <- read.csv(IN_PATH, stringsAsFactors = FALSE, colClasses = c(cond_code = "character"))

est$key    <- key(est)
limits$key <- key(limits)
if (any(!est$key %in% limits$key)) stop("극한값과 맞지 않는 조건이 있습니다 (IV 수준 확인).")

long <- est %>%
  select(cond_code, rep_id, key, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
         b3_est, b4_est) %>%
  left_join(limits %>% select(key, b3_true, b4_true, b3_truestd = b3_true_std,
                              b4_truestd = b4_true_std, b3_limit, b4_limit), by = "key") %>%
  pivot_longer(c(b3_est, b4_est, b3_true, b4_true, b3_truestd, b4_truestd, b3_limit, b4_limit),
               names_to = c("param", ".value"), names_sep = "_") %>%
  filter(is.finite(est))

summary_df <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist, param) %>%
  summarise(n          = n(),
            b_true     = first(true),
            b_true_std = first(truestd),
            b_limit    = first(limit),
            mean_est   = mean(est),
            asym_bias  = b_limit - b_true_std,
            bias_limit = mean(est - limit),
            rmse_limit = sqrt(mean((est - limit)^2)),
            bias_true  = mean(est - true),
            rmse_true  = sqrt(mean((est - true)^2)),
            .groups = "drop") %>%
  mutate(iv4_label = unname(IV4_LABEL[iv4_theta_dist])) %>%
  arrange(cond_code, param)
write.csv(summary_df, out("limit_bias_summary.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s (%d행)", out("limit_bias_summary.csv"), nrow(summary_df)))

# =============================================================================
# 3. 그래프 (09번 Bias/RMSE 그래프와 같은 배치: 행 모수, 열 IV1 × IV2, x축 IV3, 색 IV4)
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
plot_df <- summary_df %>%
  mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
         iv2 = make_num_factor(iv2_b_interval, "%.1f", decreasing = TRUE),
         iv3 = make_num_factor(iv3_b_mean,     "%.1f"),
         iv4_label = factor(iv4_label, levels = unname(IV4_LABEL[c("normal", "pos_skew",
                                                                   "neg_skew", "uniform")])),
         param = factor(param, levels = c("b3", "b4")))

common_layers <- list(
  facet_grid(rows = vars(param), cols = vars(iv1, iv2),
             labeller = labeller(iv1 = function(x) paste0("sf4=", x),
                                 iv2 = function(x) paste0("간격=", x),
                                 param = label_value)),
  scale_colour_brewer(palette = "Set1"),
  labs(x = "문항 심각도 — IV3: 경계모수 평균", colour = "능력모수 분포 (IV4)"),
  theme_bw(base_size = 10),
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
)
line_plot <- function(yvar, title, ylab, fname, zero_line = TRUE) {
  p <- ggplot(plot_df, aes(x = iv3, y = .data[[yvar]], colour = iv4_label, group = iv4_label))
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4)
  p <- p + geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
    labs(title = title, y = ylab) + common_layers
  ggsave(out(fname), p, width = 14, height = 7, dpi = 150)
  message(sprintf("저장 완료: %s", out(fname)))
}
line_plot("bias_limit", "문항 1 경계 모수 Bias (추정 - 극한값)", "Bias (극한값 기준)",
          "limit_bias_plot.png")
line_plot("rmse_limit", "문항 1 경계 모수 RMSE (극한값 기준)", "RMSE (극한값 기준)",
          "limit_rmse_plot.png", zero_line = FALSE)
line_plot("asym_bias",  "문항 1 경계 모수 극한 Bias (극한값 - 참값, 추정 척도)", "극한 Bias",
          "limit_asym_bias_plot.png")

# ── 참값·극한값·추정 평균 비교 (b3, b4 각각, 행 IV4) ──
VAL_LEVELS  <- c("참값 (추정 척도)", "극한값 (N -> inf)", "추정 평균")
val_colours <- setNames(c("grey30", "#4C9BE8", "#D1495B"), VAL_LEVELS)
for (pp in c("b3", "b4")) {
  vdf <- plot_df %>% filter(param == pp) %>%
    select(iv1, iv2, iv3, iv4_label, b_true_std, b_limit, mean_est) %>%
    pivot_longer(c(b_true_std, b_limit, mean_est), names_to = "what", values_to = "value") %>%
    mutate(what = factor(what, levels = c("b_true_std", "b_limit", "mean_est"), labels = VAL_LEVELS))
  p <- ggplot(vdf, aes(x = iv3, y = value, colour = what, linetype = what, group = what)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
    facet_grid(rows = vars(iv4_label), cols = vars(iv1, iv2),
               labeller = labeller(iv1 = function(x) paste0("sf4=", x),
                                   iv2 = function(x) paste0("간격=", x))) +
    scale_colour_manual(values = val_colours, name = NULL) +
    scale_linetype_manual(values = c("dashed", "solid", "solid"), name = NULL) +
    labs(title = sprintf("문항 1 %s: 참값 · 극한값 · 추정 평균", pp),
         x = "문항 심각도 — IV3: 경계모수 평균", y = pp) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8))
  ggsave(out(sprintf("limit_values_plot_%s.png", pp)), p, width = 14, height = 10, dpi = 150)
  message(sprintf("저장 완료: %s", out(sprintf("limit_values_plot_%s.png", pp))))
}

# =============================================================================
# 4. 정규분포 조건만 (IV4 = normal, 27 조건) — 11_analysis_normal_only.R 과 같은 폴더·배치
#   출력: output/analysis_normal/limit/
#     limit_bias_summary.csv       — 조건 × 모수(b3, b4)별 지표 (위 2단계와 같은 정의)
#     limit_bias_plot.png          — 극한값 기준 Bias  (행 모수, 열 간격 × s4, x축 IV3)
#     limit_rmse_plot.png          — 극한값 기준 RMSE
#     limit_asym_bias_plot.png     — 극한 Bias (극한값 − 참값)
#     limit_values_plot.png        — 참값·극한값·추정 평균 비교
# =============================================================================
NORM_DIR <- "output/analysis_normal/limit"
dir.create(NORM_DIR, recursive = TRUE, showWarnings = FALSE)
nout <- function(f) file.path(NORM_DIR, f)

norm_df <- summary_df %>%
  filter(iv4_theta_dist == "normal") %>%
  select(-iv4_theta_dist, -iv4_label, -b_true_std)   # 정규분포는 b_true_std = b_true
write.csv(norm_df, nout("limit_bias_summary.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s (%d행)", nout("limit_bias_summary.csv"), nrow(norm_df)))

norm_plot_df <- norm_df %>%
  mutate(f_iv1 = make_num_factor(iv1_sf4,        "%.2f"),
         f_iv2 = make_num_factor(iv2_b_interval, "%.1f", decreasing = TRUE),
         f_iv3 = make_num_factor(iv3_b_mean,     "%.1f"),
         param = factor(param, levels = c("b3", "b4")))

norm_layers <- list(
  facet_grid(rows = vars(param), cols = vars(f_iv2, f_iv1),
             labeller = labeller(f_iv2 = function(x) paste0("Interval = ", x),
                                 f_iv1 = function(x) paste0("s4 = ", x),
                                 param = label_value)),
  labs(x = "b-parameter Mean (IV3)"),
  theme_bw(base_size = 14),
  theme(strip.text = element_text(size = 11), axis.text = element_text(size = 11))
)
norm_line_plot <- function(yvar, ylab, colour, fname, zero_line = TRUE) {
  p <- ggplot(norm_plot_df, aes(x = f_iv3, y = .data[[yvar]], group = 1))
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4)
  p <- p + geom_line(linewidth = 0.7, colour = colour) + geom_point(size = 1.8, colour = colour) +
    labs(y = ylab) + norm_layers
  ggsave(nout(fname), p, width = 16, height = 7, dpi = 150)
  message(sprintf("저장 완료: %s", nout(fname)))
}
norm_line_plot("bias_limit", "Bias relative to limit (estimate - limit)", "steelblue", "limit_bias_plot.png")
norm_line_plot("rmse_limit", "RMSE relative to limit", "firebrick", "limit_rmse_plot.png", zero_line = FALSE)
norm_line_plot("asym_bias",  "Asymptotic bias (limit - true)", "darkgreen", "limit_asym_bias_plot.png")

norm_vals <- norm_plot_df %>%
  select(f_iv1, f_iv2, f_iv3, param, b_true, b_limit, mean_est) %>%
  pivot_longer(c(b_true, b_limit, mean_est), names_to = "what", values_to = "value") %>%
  mutate(what = factor(what, levels = c("b_true", "b_limit", "mean_est"),
                       labels = c("True", "Limit (N -> inf)", "Mean estimate")))
p_nv <- ggplot(norm_vals, aes(x = f_iv3, y = value, colour = what, linetype = what, group = what)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.8) +
  scale_colour_manual(values = c("True" = "grey30", "Limit (N -> inf)" = "#4C9BE8",
                                 "Mean estimate" = "#D1495B"), name = NULL) +
  scale_linetype_manual(values = c("dashed", "solid", "solid"), name = NULL) +
  labs(y = "Boundary parameter") + norm_layers +
  theme(legend.position = "bottom")
ggsave(nout("limit_values_plot.png"), p_nv, width = 16, height = 8, dpi = 150)
message(sprintf("저장 완료: %s", nout("limit_values_plot.png")))
