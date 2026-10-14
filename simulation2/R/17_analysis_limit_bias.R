# =============================================================================
# 17_analysis_limit_bias.R  (simulation2, GPCM 추정 결과)
# 채점함수 오지정 GPCM 추정에서 문항 1의 b3, b4 가 수렴하는 값(극한값)을 계산하고,
# 그 값을 기준으로 Bias / RMSE 를 계산·시각화
#
# ── 극한값(pseudo-true value) ────────────────────────────────────────────────
#   모형이 틀려도 최대우도 추정값은 참값이 아니라
#   "참 모형과 추정 모형의 KL 발산을 최소로 만드는 모수"로 수렴함 (White, 1982)
#     참 모형  : 채점함수 (0, 1, 2, s4, 4), a = 1, b = 평균 + 간격 × (−1.5, −0.5, 0.5, 1.5)
#     추정 모형: 채점함수 (0, 1, 2, 3, 4) GPCM (a, b1~b4)
#   극한값 = argmax_{a,b} ∫ φ(θ) Σ_k P_true(k|θ) log P_GPCM(k|θ; a, b) dθ ,  θ ~ N(0, 1)
#   (문항 2~20 은 올바르게 지정되어 θ 척도를 정하므로 문항 1만 따로 계산)
#
#   LIMIT_METHOD
#     "kl"         : 위 적분을 구적법으로 계산 (조건당 1초 미만, 기본값)
#                    큰 표본 MML 과 비교 시 b3·b4 차이 대부분 0.03 이내 (최대 약 0.1)
#     "mirt_large" : 조건마다 N_LARGE 명의 응답을 생성해 mirt GPCM 으로 한 번 추정한 값
#                    (실제 추정 절차와 같음, 조건당 수십 초 — 결과를 저장해 두고 재사용)
#
# ── 지표 (수렴한 반복, 조건별) ────────────────────────────────────────────────
#   bias_limit = mean(b̂ − b_limit)           — 극한값 기준 Bias (유한 표본 편향)
#   rmse_limit = sqrt(mean((b̂ − b_limit)²))  — 극한값 기준 RMSE
#   asym_bias  = b_limit − b_true             — 표본이 무한히 커도 남는 Bias (채점함수 오지정 몫)
#   bias_true / rmse_true                     — 참값 기준 (09번과 같은 값, 비교용)
#     bias_true = asym_bias + bias_limit
#
# 입력: output/analysis/rmse_bias_long.csv  (09_analysis_rmse_bias.R 출력)
# 출력: output/analysis/limit/
#   limit_values.csv         — 280 조건 전체의 극한값 (a, b1~b4) 과 참값
#   limit_bias_summary.csv   — 조건 × 모수(b3, b4)별 지표
#   limit_bias_plot.png      — 극한값 기준 Bias  (행 IV1, 열 IV2 좁음→넓음, x축 IV3, 색 b3/b4)
#   limit_rmse_plot.png      — 극한값 기준 RMSE  (같은 배치)
#   limit_asym_bias_plot.png — 극한 Bias (b_limit − b_true)  (같은 배치)
#   limit_values_plot_b3.png / _b4.png — 참값·극한값·추정 평균 비교
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

IV1_LEVELS <- c(3, 3.2, 3.33, 3.4, 3.6, 3.66, 3.8)
IV2_LEVELS <- c(1.4, 1.25, 1.2, 1, 0.8, 0.6, 0.5, 0.4, 0.25, 0.2)
IV3_LEVELS <- c(0, 0.5, 1, 1.5)
B_OFFSET   <- c(-1.5, -0.5, 0.5, 1.5)

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

QUAD   <- seq(-6, 6, length.out = 481)
QUAD_W <- dnorm(QUAD) / sum(dnorm(QUAD))

limit_kl <- function(sf4, b) {
  P   <- true_probs(QUAD, c(0, 1, 2, sf4, 4), b)
  obj <- function(par) -sum(QUAD_W * rowSums(P * log(pmax(gpcm_probs(QUAD, par), 1e-300))))
  o   <- optim(c(0, b), obj, method = "BFGS",
               control = list(reltol = 1e-14, maxit = 2000))
  if (o$convergence != 0) warning(sprintf("KL 최적화 미수렴: sf4=%g, b=%s", sf4, paste(b, collapse = ",")))
  c(a = exp(o$par[1]), b1 = o$par[2], b2 = o$par[3], b3 = o$par[4], b4 = o$par[5])
}

limit_mirt_large <- function(sf4, b, seed) {
  if (!requireNamespace("mirt", quietly = TRUE)) stop("mirt 패키지가 필요합니다.")
  set.seed(seed)
  theta <- rnorm(N_LARGE)
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

cond_grid <- expand.grid(iv3_b_mean = IV3_LEVELS, iv2_b_interval = IV2_LEVELS,
                         iv1_sf4 = IV1_LEVELS)[, 3:1]
limit_path <- out(sprintf("limit_values_%s.csv", LIMIT_METHOD))
done <- if (file.exists(limit_path)) read.csv(limit_path) else NULL
key  <- function(d) sprintf("%.4f_%.4f_%.4f", d$iv1_sf4, d$iv2_b_interval, d$iv3_b_mean)
todo <- cond_grid[!key(cond_grid) %in% if (is.null(done)) character(0) else key(done), ]

message(sprintf("극한값 계산 (%s): %d 조건 (저장된 값 %d 조건 재사용)",
                LIMIT_METHOD, nrow(todo), nrow(cond_grid) - nrow(todo)))
for (i in seq_len(nrow(todo))) {
  cg <- todo[i, ]
  b  <- cg$iv3_b_mean + cg$iv2_b_interval * B_OFFSET
  lim <- if (LIMIT_METHOD == "kl") limit_kl(cg$iv1_sf4, b) else
         limit_mirt_large(cg$iv1_sf4, b, LARGE_SEED + i)
  row <- data.frame(cg, a_true = 1,
                    b1_true = b[1], b2_true = b[2], b3_true = b[3], b4_true = b[4],
                    a_limit = lim[["a"]], b1_limit = lim[["b1"]], b2_limit = lim[["b2"]],
                    b3_limit = lim[["b3"]], b4_limit = lim[["b4"]])
  done <- bind_rows(done, row)
  if (LIMIT_METHOD == "mirt_large" || i == nrow(todo))
    write.csv(done, limit_path, row.names = FALSE)   # mirt_large 는 조건마다 저장 (중단 후 이어서)
  if (LIMIT_METHOD == "mirt_large")
    message(sprintf("  %d/%d  sf4=%g 간격=%g 평균=%g", i, nrow(todo),
                    cg$iv1_sf4, cg$iv2_b_interval, cg$iv3_b_mean))
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

match_key <- function(d) sprintf("%.4f_%.4f_%.4f", d$iv1_sf4, d$iv2_b_interval, d$iv3_b_mean)
est$key    <- match_key(est)
limits$key <- match_key(limits)
if (any(!est$key %in% limits$key)) stop("극한값과 맞지 않는 조건이 있습니다 (IV 수준 확인).")

long <- est %>%
  select(cond_code, rep_id, key, iv1_sf4, iv2_b_interval, iv3_b_mean, b3_est, b4_est) %>%
  left_join(limits %>% select(key, b3_true, b4_true, b3_limit, b4_limit), by = "key") %>%
  pivot_longer(c(b3_est, b4_est, b3_true, b4_true, b3_limit, b4_limit),
               names_to = c("param", ".value"), names_sep = "_") %>%
  filter(is.finite(est))

summary_df <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
  summarise(n          = n(),
            b_true     = first(true),
            b_limit    = first(limit),
            mean_est   = mean(est),
            asym_bias  = b_limit - b_true,
            bias_limit = mean(est - limit),
            rmse_limit = sqrt(mean((est - limit)^2)),
            bias_true  = mean(est - true),
            rmse_true  = sqrt(mean((est - true)^2)),
            .groups = "drop") %>%
  arrange(cond_code, param)
write.csv(summary_df, out("limit_bias_summary.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s (%d행)", out("limit_bias_summary.csv"), nrow(summary_df)))

# =============================================================================
# 3. 그래프 (행 IV1 × 열 IV2 좁음→넓음, x축 IV3 — 전치 그래프와 같은 배치)
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
                iv2 = make_num_factor(iv2_b_interval, "%.2f"),
                iv3 = make_num_factor(iv3_b_mean,     "%.1f"))
}
param_colours <- c(b3 = "#F2A541", b4 = "#D1495B")
common_layers <- list(
  facet_grid(rows = vars(iv1), cols = vars(iv2),
             labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                 iv2 = function(x) paste0("Int=", x))),
  scale_colour_manual(values = param_colours, name = "Boundary"),
  labs(x = "b-parameter Mean (IV3)"),
  theme_bw(base_size = 9),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))
)
plot_df <- add_factors(summary_df)

line_plot <- function(yvar, ylab, fname, zero_line = TRUE) {
  p <- ggplot(plot_df, aes(x = iv3, y = .data[[yvar]], colour = param, group = param))
  if (zero_line) p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4)
  p <- p + geom_line(linewidth = 0.7) + geom_point(size = 1.4) + labs(y = ylab) + common_layers
  ggsave(out(fname), p, width = 24, height = 18, dpi = 150)
  message(sprintf("저장 완료: %s", out(fname)))
}
line_plot("bias_limit", "Bias relative to limit (mean of estimate - limit)", "limit_bias_plot.png")
line_plot("rmse_limit", "RMSE relative to limit", "limit_rmse_plot.png", zero_line = FALSE)
line_plot("asym_bias",  "Asymptotic bias (limit - true)", "limit_asym_bias_plot.png")

# ── 참값·극한값·추정 평균 비교 (b3, b4 각각) ──
VAL_LEVELS  <- c("True", "Limit (N -> inf)", "Mean estimate")
val_colours <- setNames(c("grey30", "#4C9BE8", "#D1495B"), VAL_LEVELS)
for (pp in c("b3", "b4")) {
  vdf <- plot_df %>% filter(param == pp) %>%
    select(iv1, iv2, iv3, b_true, b_limit, mean_est) %>%
    pivot_longer(c(b_true, b_limit, mean_est), names_to = "what", values_to = "value") %>%
    mutate(what = factor(what, levels = c("b_true", "b_limit", "mean_est"), labels = VAL_LEVELS))
  p <- ggplot(vdf, aes(x = iv3, y = value, colour = what, linetype = what, group = what)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_colour_manual(values = val_colours, name = NULL) +
    scale_linetype_manual(values = c("dashed", "solid", "solid"), name = NULL) +
    labs(x = "b-parameter Mean (IV3)", y = sprintf("%s (GPCM)", pp)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out(sprintf("limit_values_plot_%s.png", pp)), p, width = 24, height = 18, dpi = 150)
  message(sprintf("저장 완료: %s", out(sprintf("limit_values_plot_%s.png", pp))))
}
