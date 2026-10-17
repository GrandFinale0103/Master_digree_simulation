# =============================================================================
# 22_plot_nrm_param_hist.R  (simulation2, NRM 추정 결과)
# 문항 1 NRM 추정값의 조건별 분포 히스토그램
#   경계모수 b4, 채점함수 s1~s3 (mirt ak1~ak3), 절편 c1~c4 (mirt d1~d4)
#   (s0 = 0, s4 = 4, c0 = 0 은 식별 제약으로 고정이라 제외)
#
# 참값 (생성 모형을 NRM 모수로 옮긴 값, a = 1)
#   s = (0, 1, 2, IV1, 4),  c_k = −Σ_{j≤k} w_j b_j  (w_j = s_j − s_{j−1}),  b4 = IV3 + 1.5 × IV2
#
# 패널 배치: 전치 비율 그래프(12_analysis_nrm_test.R 의 nrm_transposition_plot.png)와 같음
#   행 = IV1 (s4), 열 = IV2 (경계모수 간격, 왼쪽 좁음 → 오른쪽 넓음)
#   전치 그래프에서 패널 안 x축이던 IV3 (경계모수 평균) 는
#     (1) hist_overlay.png        — 한 패널에 IV3 4수준을 색으로 겹침 (선 히스토그램)
#     (2) hist_mean<IV3>.png      — IV3 수준마다 1장 (막대 히스토그램)
#
# x축과 급간: 같은 모수의 그림은 모두 같은 x축 범위와 같은 급간을 씀
#   (반복 집합 3가지, 겹친 그림, IV3 수준별 그림 4장 모두 동일 — 제외 없는 전체 반복으로 한 번 정함)
#   x축 범위: PARAM_XLIM 에 지정한 값, NA 면 참값 범위 ± 4 × (조건별 IQR 의 중앙값)
#   급간    : PARAM_BINWIDTH 에 지정한 값, NA 면 (조건별 IQR 의 중앙값) / BINS_PER_IQR
#             → 조건 하나의 분포 가운데 50% 를 약 BINS_PER_IQR 칸으로 나누는 폭
#   범위 밖 반복 수는 (2) 의 패널 글자와 param_dist_cond.csv 에 기록
#
# 표시 (IV3 수준별 그림): 검은 점선 = 참값, 빨간 실선 = 추정값 중앙값, 파란 점-선 = 추정값 평균
#   평균이 x축 범위 밖이면 가장자리에 그리고 패널 글자에 평균값을 적음 (b4 처럼 극단값이 많은 모수)
#   b4 패널 글자에는 b̂4 < 0 비율과 Δs4 ≤ 0 비율도 표기
#   겹친 그림은 IV3 4수준의 참값(색 점선)만 표시
#
# 반복 집합 (16·18·21번과 같은 기준, 반복 단위 제외)
#   all                — 아무것도 제외하지 않음 (수렴하고 b1~b4 가 계산된 반복)
#   excl_extreme       — b1~b4 중 하나라도 |b̂ − b| ≥ ERR_CUT 인 반복 제외
#   excl_extreme_sfrev — 위 + 채점함수 이상(Δs ≤ 0 또는 |Δs| < DAK_CUT) 반복 제외
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/param_hist/<반복 집합>/
#   param_dist_cond.csv                 — 조건 × 모수별 n, 참값, 평균, 분위수, 편향, 범위 밖 수
#   <모수>/hist_overlay.png, <모수>/hist_mean<IV3>.png
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
PARAMS     <- c("b4", "s1", "s2", "s3", "c1", "c2", "c3", "c4")   # 그릴 모수
BINS_PER_IQR <- 3       # 급간 자동 설정: 조건별 IQR 중앙값을 몇 칸으로 나눌지 (작을수록 급간이 넓음)
PARAM_XLIM <- list(     # 모수별 x축 범위, NA 면 자동
  b4 = c(-4, 8),
  s1 = NA, s2 = NA, s3 = NA,
  c1 = NA, c2 = NA, c3 = NA, c4 = NA
)
PARAM_BINWIDTH <- list( # 모수별 급간, NA 면 자동
  b4 = 0.5,
  s1 = NA, s2 = NA, s3 = NA,
  c1 = NA, c2 = NA, c3 = NA, c4 = NA
)
ERR_CUT    <- 3         # 극단 오차 기준 (15·16번과 같은 값)
DAK_CUT    <- 0.1       # 분모 붕괴 기준 (15·16번과 같은 값)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
BASE_DIR <- "output/nrm/param_hist"
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

PARAM_TITLE <- c(b4 = "b4 (boundary, categories 4-5)",
                 s1 = "s1 (scoring function, true = 1)",
                 s2 = "s2 (scoring function, true = 2)",
                 s3 = "s3 (scoring function, true = s4)",
                 c1 = "c1 (intercept)", c2 = "c2 (intercept)",
                 c3 = "c3 (intercept)", c4 = "c4 (intercept)")

# =============================================================================
# 1. 데이터와 참값 (16번과 같은 제외 기준)
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

bt <- sapply(B_OFFSET, function(o) rep_df$iv3_b_mean + rep_df$iv2_b_interval * o)
w  <- cbind(1, 1, rep_df$iv1_sf4 - 2, 4 - rep_df$iv1_sf4)
ct <- -t(apply(w * bt, 1, cumsum))

rep_df <- rep_df %>%
  mutate(is_extreme = pmax(abs(b1 - bt[, 1]), abs(b2 - bt[, 2]),
                           abs(b3 - bt[, 3]), abs(b4 - bt[, 4])) >= ERR_CUT,
         is_sfrev   = {
           dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
           rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
         },
         dak4 = ak4 - ak3)

# 반복 × 모수 (long)
EST_COL  <- c(b4 = "b4", s1 = "ak1", s2 = "ak2", s3 = "ak3", c1 = "d1", c2 = "d2", c3 = "d3", c4 = "d4")
true_val <- function(pp) switch(pp,
  b4 = bt[, 4], s1 = rep(1, nrow(rep_df)), s2 = rep(2, nrow(rep_df)), s3 = rep_df$iv1_sf4,
  c1 = ct[, 1], c2 = ct[, 2], c3 = ct[, 3], c4 = ct[, 4])
long <- bind_rows(lapply(PARAMS, function(pp) {
  data.frame(cond_code = rep_df$cond_code, iv1_sf4 = rep_df$iv1_sf4,
             iv2_b_interval = rep_df$iv2_b_interval, iv3_b_mean = rep_df$iv3_b_mean,
             param = pp, est = rep_df[[EST_COL[[pp]]]], true = true_val(pp),
             dak4 = rep_df$dak4, is_extreme = rep_df$is_extreme, is_sfrev = rep_df$is_sfrev,
             stringsAsFactors = FALSE)
})) %>% filter(!is.nan(est))

SETS <- list(
  all                = function(d) rep(TRUE, nrow(d)),
  excl_extreme       = function(d) !d$is_extreme,
  excl_extreme_sfrev = function(d) !d$is_extreme & !d$is_sfrev
)
SET_DESC <- c(all = "no exclusion",
              excl_extreme = sprintf("excluding extreme errors (|error| >= %g)", ERR_CUT),
              excl_extreme_sfrev = sprintf("excluding extreme errors and SF reversals / near-zero ds (|ds| < %g)", DAK_CUT))

# =============================================================================
# 2. 그래프 도우미
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
                iv2 = make_num_factor(iv2_b_interval, "%.2f"),   # 좁음 → 넓음 (전치 그래프와 같음)
                iv3 = make_num_factor(iv3_b_mean,     "%.1f"))
}
grid_facet <- facet_grid(rows = vars(iv1), cols = vars(iv2), scales = "free_y",
                         labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                             iv2 = function(x) paste0("Int=", x)))
base_theme <- list(
  theme_bw(base_size = 9),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7))
)
IV3_COLOURS <- c("0.0" = "#4C9BE8", "0.5" = "#59A14F", "1.0" = "#F2A541", "1.5" = "#D1495B")

# 모수별 x축 범위와 급간: 제외 없는 전체 반복으로 한 번만 정해 모든 그림에 같이 씀
AXIS <- lapply(setNames(PARAMS, PARAMS), function(pp) {
  d   <- long %>% filter(param == pp, is.finite(est))
  iqr <- d %>% group_by(cond_code) %>% summarise(i = IQR(est), .groups = "drop") %>% pull(i) %>% median()
  r   <- range(d$true)
  xl  <- if (!anyNA(PARAM_XLIM[[pp]])) PARAM_XLIM[[pp]] else c(r[1] - 4 * iqr, r[2] + 4 * iqr)
  bw  <- if (!anyNA(PARAM_BINWIDTH[[pp]])) PARAM_BINWIDTH[[pp]] else iqr / BINS_PER_IQR
  list(xlim = xl, bw = bw)
})

# =============================================================================
# 3. 반복 집합 × 모수별 표와 그래프
# =============================================================================
for (st in names(SETS)) {
  set_dir <- file.path(BASE_DIR, st)
  dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
  ds <- add_factors(long[SETS[[st]](long), ])

  dist_all <- list()
  for (pp in PARAMS) {
    d    <- ds %>% filter(param == pp)
    xlim <- AXIS[[pp]]$xlim
    bw   <- AXIS[[pp]]$bw
    out  <- function(f) file.path(set_dir, pp, f)
    dir.create(file.path(set_dir, pp), recursive = TRUE, showWarnings = FALSE)

    # ── 조건별 분포 요약 ──
    dist_cond <- d %>%
      group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
      summarise(n          = n(),
                true       = first(true),
                mean_est   = mean(est[is.finite(est)]),
                q05 = quantile(est, 0.05), q25 = quantile(est, 0.25), median_est = median(est),
                q75 = quantile(est, 0.75), q95 = quantile(est, 0.95),
                share_below_true = mean(est < true),
                share_negative   = mean(est < 0),
                share_ds4_le0    = mean(dak4 <= 0),
                n_outside_xlim   = sum(est < xlim[1] | est > xlim[2]),
                .groups = "drop") %>%
      mutate(mean_bias = mean_est - true, median_bias = median_est - true,
             xlim_low = xlim[1], xlim_high = xlim[2], binwidth = bw)
    dist_all[[pp]] <- dist_cond

    inside <- d %>% filter(est >= xlim[1], est <= xlim[2])
    lines  <- add_factors(dist_cond)

    # ── (1) IV3 4수준을 한 패널에 겹침 ──
    p_ov <- ggplot(inside, aes(x = est, colour = iv3)) +
      geom_freqpoly(binwidth = bw, boundary = xlim[1], linewidth = 0.6) +
      geom_vline(data = lines, aes(xintercept = true, colour = iv3),
                 linetype = "dashed", linewidth = 0.4) +
      grid_facet +
      scale_colour_manual(values = IV3_COLOURS, name = "b-parameter Mean (IV3)") +
      coord_cartesian(xlim = xlim) +
      labs(title = sprintf("Item 1 NRM %s: distribution of estimates (dashed = true) - %s",
                           PARAM_TITLE[[pp]], SET_DESC[[st]]),
           x = sprintf("%s estimate (bin width %.2f)", pp, bw), y = "Count") +
      base_theme
    ggsave(out("hist_overlay.png"), p_ov, width = 24, height = 18, dpi = 150)

    # ── (2) IV3 수준마다 1장 ──
    for (m in sort(unique(d$iv3_b_mean))) {
      h  <- inside %>% filter(iv3_b_mean == m)
      vl <- lines %>% filter(iv3_b_mean == m)
      clamp <- function(x) pmin(pmax(x, xlim[1]), xlim[2])
      vl_long <- bind_rows(vl %>% transmute(iv1, iv2, value = true,              what = "True value"),
                           vl %>% transmute(iv1, iv2, value = clamp(median_est), what = "Median of estimates"),
                           vl %>% transmute(iv1, iv2, value = clamp(mean_est),   what = "Mean of estimates")) %>%
        mutate(what = factor(what, levels = c("True value", "Median of estimates", "Mean of estimates")))
      lab <- vl %>% mutate(
        mean_txt = ifelse(mean_est < xlim[1] | mean_est > xlim[2],
                          sprintf("mean %.2f (off-scale)", mean_est), sprintf("mean %.2f", mean_est)),
        label = if (pp == "b4")
          sprintf("n=%d, outside: %d\n%s, median %.2f\nb4<0: %.0f%%, ds4<=0: %.0f%%",
                  n, n_outside_xlim, mean_txt, median_est, 100 * share_negative, 100 * share_ds4_le0) else
          sprintf("n=%d, outside: %d\n%s, median %.2f\nbelow true: %.0f%%",
                  n, n_outside_xlim, mean_txt, median_est, 100 * share_below_true))
      p_h <- ggplot(h, aes(x = est)) +
        geom_histogram(binwidth = bw, boundary = xlim[1], fill = "grey70", colour = "grey40", linewidth = 0.2) +
        geom_vline(data = vl_long, aes(xintercept = value, colour = what, linetype = what), linewidth = 0.6) +
        geom_text(data = lab, aes(x = xlim[1], y = Inf, label = label),
                  hjust = 0, vjust = 1.2, size = 2.2, lineheight = 0.9, inherit.aes = FALSE) +
        grid_facet +
        scale_colour_manual(values = c("True value" = "black", "Median of estimates" = "firebrick",
                                       "Mean of estimates" = "#2C6FB7"), name = NULL) +
        scale_linetype_manual(values = c("True value" = "dashed", "Median of estimates" = "solid",
                                         "Mean of estimates" = "dotdash"), name = NULL) +
        coord_cartesian(xlim = xlim) +
        labs(title = sprintf("Item 1 NRM %s, b-mean (IV3) = %.1f - %s", PARAM_TITLE[[pp]], m, SET_DESC[[st]]),
             x = sprintf("%s estimate (bin width %.2f)", pp, bw), y = "Count") +
        base_theme
      ggsave(out(sprintf("hist_mean%.1f.png", m)), p_h, width = 24, height = 18, dpi = 150)
    }
  }
  write.csv(bind_rows(dist_all), file.path(set_dir, "param_dist_cond.csv"), row.names = FALSE)
  message(sprintf("[%s] 반복 %d개, 모수 %d개 → %s/", st,
                  sum(ds$param == PARAMS[1]), length(PARAMS), set_dir))
}
