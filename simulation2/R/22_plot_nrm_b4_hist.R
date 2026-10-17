# =============================================================================
# 22_plot_nrm_b4_hist.R  (simulation2, NRM 추정 결과)
# 문항 1 경계모수 b4 추정값(b̂4)의 분포를 조건별 히스토그램으로 그림
#
# 패널 배치: 전치 비율 그래프(12_analysis_nrm_test.R 의 nrm_transposition_plot.png)와 같음
#   행 = IV1 (s4), 열 = IV2 (경계모수 간격, 왼쪽 좁음 → 오른쪽 넓음)
#   전치 그래프에서 패널 안 x축이던 IV3 (경계모수 평균) 는
#     (1) b4_hist_overlay.png        — 한 패널에 IV3 4수준을 색으로 겹침 (선 히스토그램)
#     (2) b4_hist_mean<IV3>.png      — IV3 수준마다 1장 (막대 히스토그램)
#   로 나타냄
#
# 표시
#   점선 세로선 = 참값 b4 (= IV3 + 1.5 × IV2),  실선 세로선 = b̂4 의 중앙값
#   x축은 XLIM 으로 자름 (b̂4 는 분모 α(s4 − s3) 가 0 근처면 매우 크거나 음수로 튐)
#   → 범위 밖 반복 수는 (2) 의 패널 글자와 b4_dist_cond.csv 에 기록
#
# 반복 집합 (16·18·21번과 같은 기준, 반복 단위 제외)
#   all                — 아무것도 제외하지 않음 (수렴하고 b1~b4 가 계산된 반복)
#   excl_extreme       — b1~b4 중 하나라도 |b̂ − b| ≥ ERR_CUT 인 반복 제외
#   excl_extreme_sfrev — 위 + 채점함수 이상(Δs ≤ 0 또는 |Δs| < DAK_CUT) 반복 제외
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/b4_hist/<반복 집합>/
#   b4_hist_overlay.png, b4_hist_mean<IV3>.png
#   b4_dist_cond.csv — 조건별 n, 참값, 평균, 분위수(5·25·50·75·95%), 음수 비율, 범위 밖 수, Δs4 ≤ 0 비율
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
XLIM     <- c(-4, 8)   # 히스토그램 x축 범위
BINWIDTH <- 0.25       # 막대 폭
ERR_CUT  <- 3          # 극단 오차 기준 (15·16번과 같은 값)
DAK_CUT  <- 0.1        # 분모 붕괴 기준 (15·16번과 같은 값)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
BASE_DIR <- "output/nrm/b4_hist"
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

# =============================================================================
# 1. 데이터 (16번과 같은 제외 기준)
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

bt <- sapply(B_OFFSET, function(o) rep_df$iv3_b_mean + rep_df$iv2_b_interval * o)
rep_df <- rep_df %>%
  mutate(b4_true    = bt[, 4],
         dak4       = ak4 - ak3,
         is_extreme = pmax(abs(b1 - bt[, 1]), abs(b2 - bt[, 2]),
                           abs(b3 - bt[, 3]), abs(b4 - bt[, 4])) >= ERR_CUT,
         is_sfrev   = {
           dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
           rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
         }) %>%
  filter(!is.nan(b4))

SETS <- list(
  all                = rep(TRUE, nrow(rep_df)),
  excl_extreme       = !rep_df$is_extreme,
  excl_extreme_sfrev = !rep_df$is_extreme & !rep_df$is_sfrev
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

# =============================================================================
# 3. 반복 집합별 표와 그래프
# =============================================================================
for (st in names(SETS)) {
  out_dir <- file.path(BASE_DIR, st)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)
  d <- add_factors(rep_df[SETS[[st]], ])

  # ── 조건별 분포 요약 ──
  dist_cond <- d %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    summarise(n          = n(),
              b4_true    = first(b4_true),
              mean_b4    = mean(b4[is.finite(b4)]),
              q05 = quantile(b4, 0.05), q25 = quantile(b4, 0.25), median_b4 = median(b4),
              q75 = quantile(b4, 0.75), q95 = quantile(b4, 0.95),
              share_negative  = mean(b4 < 0),
              share_below_true = mean(b4 < b4_true),
              n_outside_xlim  = sum(b4 < XLIM[1] | b4 > XLIM[2]),
              share_ds4_le0   = mean(dak4 <= 0),
              .groups = "drop") %>%
    mutate(median_bias = median_b4 - b4_true, mean_bias = mean_b4 - b4_true)
  write.csv(dist_cond, out("b4_dist_cond.csv"), row.names = FALSE)

  inside <- d %>% filter(b4 >= XLIM[1], b4 <= XLIM[2])
  lines  <- add_factors(dist_cond)

  # ── (1) IV3 4수준을 한 패널에 겹침 ──
  p_ov <- ggplot(inside, aes(x = b4, colour = iv3)) +
    geom_freqpoly(binwidth = BINWIDTH, boundary = 0, linewidth = 0.5) +
    geom_vline(data = lines, aes(xintercept = b4_true, colour = iv3),
               linetype = "dashed", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    grid_facet +
    scale_colour_manual(values = IV3_COLOURS, name = "b-parameter Mean (IV3)") +
    coord_cartesian(xlim = XLIM) +
    labs(title = sprintf("Item 1 NRM: distribution of b4 estimates by condition (dashed = true b4) - %s", SET_DESC[[st]]),
         x = "b4 estimate", y = "Count") +
    base_theme
  ggsave(out("b4_hist_overlay.png"), p_ov, width = 24, height = 18, dpi = 150)

  # ── (2) IV3 수준마다 1장 ──
  for (m in sort(unique(d$iv3_b_mean))) {
    h  <- inside %>% filter(iv3_b_mean == m)
    vl <- lines %>% filter(iv3_b_mean == m)
    vl_long <- bind_rows(vl %>% transmute(iv1, iv2, value = b4_true,  what = "True b4"),
                         vl %>% transmute(iv1, iv2, value = pmin(pmax(median_b4, XLIM[1]), XLIM[2]),
                                          what = "Median of estimates"))
    lab <- vl %>% mutate(label = sprintf("n=%d, outside: %d\nb4<0: %.0f%%, ds4<=0: %.0f%%",
                                         n, n_outside_xlim, 100 * share_negative, 100 * share_ds4_le0))
    p_h <- ggplot(h, aes(x = b4)) +
      geom_histogram(binwidth = BINWIDTH, boundary = 0, fill = "grey70", colour = "grey40", linewidth = 0.2) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
      geom_vline(data = vl_long, aes(xintercept = value, colour = what, linetype = what), linewidth = 0.6) +
      geom_text(data = lab, aes(x = XLIM[1], y = Inf, label = label),
                hjust = 0, vjust = 1.2, size = 2.2, lineheight = 0.9, inherit.aes = FALSE) +
      grid_facet +
      scale_colour_manual(values = c("True b4" = "black", "Median of estimates" = "firebrick"), name = NULL) +
      scale_linetype_manual(values = c("True b4" = "dashed", "Median of estimates" = "solid"), name = NULL) +
      coord_cartesian(xlim = XLIM) +
      labs(title = sprintf("Item 1 NRM: distribution of b4 estimates, b-mean (IV3) = %.1f - %s", m, SET_DESC[[st]]),
           x = "b4 estimate", y = "Count") +
      base_theme
    ggsave(out(sprintf("b4_hist_mean%.1f.png", m)), p_h, width = 24, height = 18, dpi = 150)
  }
  message(sprintf("[%s] 반복 %d개 → %s/", st, nrow(d), out_dir))
}
