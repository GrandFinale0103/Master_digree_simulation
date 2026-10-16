# =============================================================================
# 21_analysis_nrm_median_bias.R  (simulation2, NRM 추정 결과)
# 문항 1 경계모수 b1~b4 의 오차(추정값 − 참값)를 중앙값으로 요약해 편향 여부를 판단하고,
# 오차 분포를 히스토그램으로 그림
#
# ── 중앙값 편향 ──────────────────────────────────────────────────────────────
#   Bias          = mean(b̂ − b)       — 평균 (극단값에 크게 흔들림)
#   median_bias   = median(b̂ − b)     — 중앙값 (극단값에 거의 영향 없음)
#   95% 신뢰구간  : 순서통계량 기반 (분포 가정 없음, 이항분포로 순위를 정함)
#   편향 검정     : 부호 검정 H0 median = 0  ⇔  P(b̂ − b > 0) = 0.5  (이항검정, 양측)
#                   조건이 많으므로 BH(FDR) 보정, 보정 p < .05 이면 편향 있음으로 판정
#   판정          : over (중앙값 > 0, 유의) / under (중앙값 < 0, 유의) / no bias (유의하지 않음)
#
# ── 반복 집합 (16번과 같은 기준, 반복 단위 제외) ─────────────────────────────
#   all                — 아무것도 제외하지 않음 (수렴하고 b1~b4 가 계산된 반복 전체)
#   excl_extreme       — b1~b4 중 하나라도 |b̂ − b| ≥ ERR_CUT 인 반복 제외
#   excl_extreme_sfrev — 위 + 채점함수 이상(Δs ≤ 0 또는 |Δs| < DAK_CUT) 반복 제외
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/median_bias/
#   median_bias_cond.csv         — 집합 × 조건 × 경계모수: n, 평균 Bias, 중앙값, 95% CI, 양수 비율,
#                                  부호 검정 p, BH 보정 p, 판정, MdAE, RMSE
#   median_bias_overall.csv      — 집합 × 경계모수: 판정별 조건 수와 전체 오차 요약
#   median_bias_summary.txt      — 위 두 표의 요약
#   <집합>/median_bias_plot.png  — 조건별 중앙값 편향 + 95% CI (색 = 판정), 경계모수마다 패널 행
#   <집합>/mean_vs_median_plot.png — 평균 Bias 와 중앙값 편향 비교 (경계모수마다 1행)
#   <집합>/hist_b<k>_mean<IV3>.png — 조건별 오차 히스토그램 (행 IV1 × 열 IV2, IV3 수준마다 1장)
#                                    세로선: 0(검정), 평균(파랑), 중앙값(빨강)
#                                    x축은 ±HIST_XLIM 로 자르고, 범위 밖 반복 수를 패널에 표기
#   <집합>/hist_pooled.png       — 경계모수별 전체 조건 합친 오차 히스토그램 (pseudo-log 축, 자르지 않음)
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
ERR_CUT   <- 3     # 극단 오차: |b̂ − b| ≥ ERR_CUT   (15·16번과 같은 값)
DAK_CUT   <- 0.1   # 분모 붕괴: |Δs_k| < DAK_CUT      (15·16번과 같은 값)
ALPHA     <- 0.05  # 부호 검정 유의수준 (BH 보정 후)
HIST_XLIM <- 6     # 조건별 히스토그램 x축 범위 ±HIST_XLIM
HIST_BINS <- 40
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
OUT_DIR <- "output/nrm/median_bias"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)
B_PARAMS <- paste0("b", 1:4)
JUDGE_LEVELS  <- c("over", "under", "no bias")
JUDGE_COLOURS <- c(over = "#D1495B", under = "#4C9BE8", `no bias` = "grey55")

# =============================================================================
# 1. 데이터: 반복 × 경계모수
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

bt  <- sapply(B_OFFSET, function(o) rep_df$iv3_b_mean + rep_df$iv2_b_interval * o)
err <- as.matrix(rep_df[, B_PARAMS]) - bt
dak <- cbind(rep_df$ak1 - rep_df$ak0, rep_df$ak2 - rep_df$ak1,
             rep_df$ak3 - rep_df$ak2, rep_df$ak4 - rep_df$ak3)
rep_df$is_extreme <- apply(abs(err) >= ERR_CUT, 1, any)
rep_df$is_sfrev   <- rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0

long <- bind_rows(lapply(1:4, function(k) {
  data.frame(cond_code = rep_df$cond_code, rep_id = rep_df$rep_id,
             iv1_sf4 = rep_df$iv1_sf4, iv2_b_interval = rep_df$iv2_b_interval,
             iv3_b_mean = rep_df$iv3_b_mean, param = B_PARAMS[k], err = err[, k],
             is_extreme = rep_df$is_extreme, is_sfrev = rep_df$is_sfrev,
             stringsAsFactors = FALSE)
})) %>% filter(!is.nan(err))          # NaN(0/0) 제외, ±Inf 는 중앙값·부호에는 그대로 사용

SETS <- list(
  all                = function(d) rep(TRUE, nrow(d)),
  excl_extreme       = function(d) !d$is_extreme,
  excl_extreme_sfrev = function(d) !d$is_extreme & !d$is_sfrev
)
SET_DESC <- c(all = "no exclusion",
              excl_extreme = sprintf("excluding extreme errors (|error| >= %g)", ERR_CUT),
              excl_extreme_sfrev = sprintf("excluding extreme errors and SF reversals / near-zero ds (|ds| < %g)", DAK_CUT))

# =============================================================================
# 2. 중앙값 편향 (조건 × 경계모수)
# =============================================================================
# 중앙값의 순서통계량 신뢰구간: 순위 l, u 는 이항분포 Bin(n, 0.5) 분위수
median_ci <- function(x, level = 0.95) {
  x <- sort(x); n <- length(x)
  if (n < 6) return(c(NA_real_, NA_real_))
  l <- qbinom((1 - level) / 2, n, 0.5)
  u <- n - l + 1
  c(x[max(l, 1)], x[min(u, n)])
}
sign_p <- function(x) {
  x <- x[x != 0]
  if (length(x) == 0) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5)$p.value
}

cond_df <- bind_rows(lapply(names(SETS), function(st) {
  long[SETS[[st]](long), ] %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
    summarise(n           = n(),
              bias_mean   = mean(err[is.finite(err)]),
              median_bias = median(err),
              ci_low      = median_ci(err)[1],
              ci_high     = median_ci(err)[2],
              prop_pos    = mean(err > 0),
              p_sign      = sign_p(err),
              mdae        = median(abs(err)),
              rmse        = sqrt(mean(err[is.finite(err)]^2)),
              .groups = "drop") %>%
    mutate(set = st, .before = 1)
})) %>%
  group_by(set) %>%
  mutate(p_sign_BH = p.adjust(p_sign, method = "BH")) %>%
  ungroup() %>%
  mutate(judgement = case_when(is.na(p_sign_BH) | p_sign_BH >= ALPHA ~ "no bias",
                               median_bias > 0                       ~ "over",
                               TRUE                                  ~ "under"),
         judgement = factor(judgement, levels = JUDGE_LEVELS),
         set = factor(set, levels = names(SETS))) %>%
  arrange(set, cond_code, param)
write.csv(cond_df, file.path(OUT_DIR, "median_bias_cond.csv"), row.names = FALSE)

overall_df <- cond_df %>%
  group_by(set, param) %>%
  summarise(n_conditions   = n(),
            n_over         = sum(judgement == "over"),
            n_under        = sum(judgement == "under"),
            n_no_bias      = sum(judgement == "no bias"),
            mean_of_median_bias = mean(median_bias),
            mean_of_mean_bias   = mean(bias_mean),
            .groups = "drop")
write.csv(overall_df, file.path(OUT_DIR, "median_bias_overall.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "median_bias_summary.txt"))
tryCatch({
  cat("NRM 문항 1 경계모수 오차(추정 - 참값)의 중앙값 편향\n")
  cat(sprintf("편향 판정: 부호 검정(H0 median = 0), 반복 집합별 BH 보정 p < %.2f\n\n", ALPHA))
  cat("반복 집합 × 경계모수별 판정 조건 수와 평균·중앙값 편향의 조건 평균\n")
  print(as.data.frame(overall_df), digits = 4, row.names = FALSE)
  cat("\n반복 집합별 남은 반복 수\n")
  for (st in names(SETS)) cat(sprintf("  %-20s %d\n", st, sum(SETS[[st]](rep_df))))
}, finally = sink())
message(sprintf("저장 완료: %s/median_bias_cond.csv, median_bias_overall.csv, median_bias_summary.txt", OUT_DIR))

# =============================================================================
# 3. 그래프
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
                iv2 = make_num_factor(iv2_b_interval, "%.2f"),   # 좁음 → 넓음
                iv3 = make_num_factor(iv3_b_mean,     "%.1f"))
}
base_theme <- list(
  theme_bw(base_size = 9),
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7), axis.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))
)

for (st in names(SETS)) {
  set_dir <- file.path(OUT_DIR, st)
  dir.create(set_dir, recursive = TRUE, showWarnings = FALSE)
  sd_ <- add_factors(cond_df %>% filter(set == st))

  # ── 조건별 중앙값 편향 + 95% CI (행: 경계모수 × IV1, 열: IV2) ──
  p_med <- ggplot(sd_, aes(x = iv3, y = median_bias, colour = judgement)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.25, linewidth = 0.4) +
    geom_point(size = 1.4) +
    facet_grid(rows = vars(param, iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_colour_manual(values = JUDGE_COLOURS, drop = FALSE,
                        name = sprintf("Sign test (BH p < %.2f)", ALPHA)) +
    labs(title = sprintf("Item 1 NRM boundary: median of (estimate - true) with 95%% CI - %s", SET_DESC[[st]]),
         x = "b-parameter Mean (IV3)", y = "Median bias") +
    base_theme
  ggsave(file.path(set_dir, "median_bias_plot.png"), p_med, width = 24, height = 40, dpi = 120,
         limitsize = FALSE)

  # ── 평균 Bias vs 중앙값 편향 ──
  mm <- sd_ %>% select(iv1, iv2, iv3, param, bias_mean, median_bias) %>%
    pivot_longer(c(bias_mean, median_bias), names_to = "stat", values_to = "value") %>%
    mutate(stat = factor(stat, levels = c("bias_mean", "median_bias"),
                         labels = c("Mean bias", "Median bias")))
  p_mm <- ggplot(mm, aes(x = iv3, y = value, colour = stat, group = stat)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.2) +
    facet_grid(rows = vars(param, iv1), cols = vars(iv2), scales = "free_y",
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_colour_manual(values = c("Mean bias" = "steelblue", "Median bias" = "firebrick"), name = NULL) +
    labs(title = sprintf("Item 1 NRM boundary: mean vs median bias - %s", SET_DESC[[st]]),
         x = "b-parameter Mean (IV3)", y = "Bias") +
    base_theme
  ggsave(file.path(set_dir, "mean_vs_median_plot.png"), p_mm, width = 24, height = 40, dpi = 120,
         limitsize = FALSE)

  # ── 조건별 오차 히스토그램 (경계모수 × IV3 수준마다 1장) ──
  hd <- add_factors(long[SETS[[st]](long), ])
  stat_df <- add_factors(cond_df %>% filter(set == st)) %>%
    select(iv1, iv2, iv3, iv3_b_mean, param, bias_mean, median_bias)
  out_df <- hd %>% group_by(iv1, iv2, iv3_b_mean, param) %>%
    summarise(n_out = sum(abs(err) > HIST_XLIM), n = n(), .groups = "drop")
  for (pp in B_PARAMS) for (m in sort(unique(hd$iv3_b_mean))) {
    h  <- hd %>% filter(param == pp, iv3_b_mean == m, abs(err) <= HIST_XLIM)
    sl <- stat_df %>% filter(param == pp, iv3_b_mean == m) %>%
      pivot_longer(c(bias_mean, median_bias), names_to = "stat", values_to = "value") %>%
      mutate(stat = factor(stat, levels = c("bias_mean", "median_bias"),
                           labels = c("Mean", "Median")),
             value = pmax(pmin(value, HIST_XLIM), -HIST_XLIM))   # 범위 밖 평균은 가장자리에 표시
    ol <- out_df %>% filter(param == pp, iv3_b_mean == m) %>%
      mutate(label = sprintf("n=%d, outside +/-%g: %d", n, HIST_XLIM, n_out))
    p_h <- ggplot(h, aes(x = err)) +
      geom_histogram(bins = HIST_BINS, boundary = 0, fill = "grey70", colour = "grey40", linewidth = 0.2) +
      geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
      geom_vline(data = sl, aes(xintercept = value, colour = stat), linewidth = 0.6) +
      geom_text(data = ol, aes(x = -HIST_XLIM, y = Inf, label = label),
                hjust = 0, vjust = 1.5, size = 2.3, inherit.aes = FALSE) +
      facet_grid(rows = vars(iv1), cols = vars(iv2), scales = "free_y",
                 labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                     iv2 = function(x) paste0("Int=", x))) +
      scale_colour_manual(values = c(Mean = "steelblue", Median = "firebrick"), name = NULL) +
      coord_cartesian(xlim = c(-HIST_XLIM, HIST_XLIM)) +
      labs(title = sprintf("Item 1 NRM %s: distribution of (estimate - true), b-mean = %.1f - %s",
                           pp, m, SET_DESC[[st]]),
           x = "estimate - true", y = "Count") +
      base_theme
    ggsave(file.path(set_dir, sprintf("hist_%s_mean%.1f.png", pp, m)), p_h,
           width = 24, height = 18, dpi = 150)
  }

  # ── 전체 조건 합친 히스토그램 (pseudo-log 축, 자르지 않음) ──
  hp <- hd %>% filter(is.finite(err))
  p_pool <- ggplot(hp, aes(x = err)) +
    geom_histogram(bins = 100, fill = "grey70", colour = "grey40", linewidth = 0.2) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
    facet_wrap(~ param, ncol = 2, scales = "free_y") +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 0.5),
                       breaks = c(-1000, -100, -10, -3, -1, 0, 1, 3, 10, 100, 1000)) +
    labs(title = sprintf("Item 1 NRM boundary: (estimate - true), all conditions pooled - %s", SET_DESC[[st]]),
         x = "estimate - true (pseudo-log axis)", y = "Count") +
    theme_bw(base_size = 12)
  ggsave(file.path(set_dir, "hist_pooled.png"), p_pool, width = 12, height = 8, dpi = 150)

  message(sprintf("[%s] 그래프 저장 → %s/", st, set_dir))
}
