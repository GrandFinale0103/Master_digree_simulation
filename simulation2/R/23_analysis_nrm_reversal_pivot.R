# =============================================================================
# 23_analysis_nrm_reversal_pivot.R  (simulation2, NRM 추정 결과)
# 채점함수가 전치될 때 절편모수와 경계모수가 어떻게 추정되는지 요약 (문항 1)
#
# 경계모수 k (범주 k 와 k+1 사이, 1부터)
#   분모 D_k = α (s_k − s_{k−1})        — 기울기 × 채점함수 간격   (mirt: a1 · (ak_k − ak_{k−1}))
#   분자 N_k = c_{k−1} − c_k            — 절편 간격의 부호 반대      (mirt: d_{k−1} − d_k)
#   b_k = N_k / D_k
#   인접 범주 로짓: logit(범주 k+1 vs k | θ) = D_k θ − N_k
#
# ── 1) 채점함수 상태별 요약 (reversal_table_cond.csv) ─────────────────────────
#   조건 × 경계 × 상태(서열 유지 D_k > 0 / 전치 D_k ≤ 0) 마다
#     n, 절편 간격 N_k 의 중앙값과 음수 비율, b̂_k 의 중앙값·음수 비율·참값보다 작은 비율
#
# ── 2) 채점함수와 절편의 상관 (correlation_cond.csv) ───────────────────────────
#   조건 × 경계마다 반복들에 걸친 상관
#     cor_N_D      : 분자 N_k 와 분모 D_k (Pearson / Spearman)
#     cor_ds_dc    : 채점함수 간격 Δs_k = s_k − s_{k−1} 과 절편 간격 Δc_k = c_k − c_{k−1}
#     cor_s_c      : 채점함수 값 s_k 와 절편 값 c_k   (k = 1..3; k = 4 는 s4 고정이라 NA)
#     cor_slo_chi  : 아래쪽 채점함수 s_{k−1} 와 위쪽 절편 c_k     (b4 면 s3 와 c4)
#     cor_slo_clo  : 아래쪽 채점함수 s_{k−1} 와 아래쪽 절편 c_{k−1} (b4 면 s3 와 c3)
#   피벗 분석: N_k = a + θ* · D_k 회귀
#     θ*        : 직선이 돌아가는 축 (데이터가 인접 범주 비율을 정해 주는 θ)
#     L*        : θ* 에서의 인접 범주 로짓 = D_k θ* − N_k  (중앙값, IQR)
#     → b_k = N_k / D_k ≈ θ* − L* / D_k  : D_k ≤ 0 이면 b_k 는 항상 θ* 아래로 떨어짐
#
# 반복 집합: 수렴하고 b1~b4 가 계산된 반복 전체 (제외 없음 — 채점함수 이상 반복을 빼면
#           전치 반복 자체가 사라지므로)
# 상관과 회귀는 유한한 값만 사용, Pearson 은 극단값에 민감하므로 Spearman 도 함께 기록
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/reversal_pivot/
#   reversal_table_cond.csv  — 1) 조건 × 경계 × 채점함수 상태
#   correlation_cond.csv     — 2) 조건 × 경계 상관·피벗
#   reversal_summary.txt     — 경계별 전체 합계와 s4 × 경계모수 평균(IV3) 별 요약
#   corr_N_D_<경계>.png      — 조건별 cor(N_k, D_k) (전치 그래프 배치: 행 IV1, 열 IV2, x축 IV3)
#   reversed_share_b_neg_<경계>.png — 조건별 전치 비율과, 전치 반복 중 b̂ < 0 비율
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
BOUNDARIES <- 4          # 요약할 경계모수 (예: 1:4 로 바꾸면 b1~b4 모두)
MIN_N_CORR <- 10         # 상관·회귀를 계산할 최소 반복 수
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
OUT_DIR <- "output/nrm/reversal_pivot"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

# =============================================================================
# 1. 데이터: 반복 × 경계 (long)
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

long <- bind_rows(lapply(BOUNDARIES, function(k) {
  s_hi <- rep_df[[paste0("ak", k)]]; s_lo <- rep_df[[paste0("ak", k - 1)]]
  c_hi <- rep_df[[paste0("d",  k)]]; c_lo <- rep_df[[paste0("d",  k - 1)]]
  data.frame(cond_code = rep_df$cond_code, rep_id = rep_df$rep_id,
             iv1_sf4 = rep_df$iv1_sf4, iv2_b_interval = rep_df$iv2_b_interval,
             iv3_b_mean = rep_df$iv3_b_mean,
             boundary = paste0("b", k),
             b_true = rep_df$iv3_b_mean + rep_df$iv2_b_interval * B_OFFSET[k],
             s_k = if (k <= 3) s_hi else NA_real_,      # s4 는 식별 제약으로 고정
             s_lo = if (k >= 2) s_lo else NA_real_,     # s0 는 식별 제약으로 고정
             c_k = c_hi, c_lo = c_lo,
             ds  = s_hi - s_lo,                         # 채점함수 간격 Δs_k
             dc  = c_hi - c_lo,                         # 절편 간격 Δc_k
             D   = rep_df$a1 * (s_hi - s_lo),           # 분모
             N   = c_lo - c_hi,                         # 분자
             b_hat = rep_df[[paste0("b", k)]],
             stringsAsFactors = FALSE)
})) %>%
  mutate(sf_state = ifelse(D > 0, "ordered (D > 0)", "reversed (D <= 0)"))

# =============================================================================
# 2. 채점함수 상태별 요약표
# =============================================================================
state_summary <- function(df) df %>%
  summarise(n                  = n(),
            median_N           = median(N),
            share_N_negative   = mean(N < 0),
            median_D           = median(D),
            median_b           = median(b_hat),
            share_b_negative   = mean(b_hat < 0),
            share_b_below_true = mean(b_hat < b_true),
            median_b_bias      = median(b_hat - b_true),
            .groups = "drop")

rev_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, b_true, sf_state) %>%
  state_summary() %>%
  group_by(cond_code, boundary) %>%
  mutate(share_of_condition = n / sum(n)) %>%
  ungroup() %>%
  arrange(cond_code, boundary, sf_state)
write.csv(rev_cond, out("reversal_table_cond.csv"), row.names = FALSE)

# =============================================================================
# 3. 상관과 피벗
# =============================================================================
safe_cor <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < MIN_N_CORR || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok], method = method)
}
pivot_fit <- function(N, D) {
  ok <- is.finite(N) & is.finite(D)
  if (sum(ok) < MIN_N_CORR) return(c(NA, NA, NA, NA))
  f  <- lm(N[ok] ~ D[ok]); th <- unname(coef(f)[2])
  L  <- D[ok] * th - N[ok]
  c(unname(coef(f)[1]), th, median(L), IQR(L))
}

corr_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, b_true) %>%
  summarise(n                 = n(),
            share_reversed    = mean(D <= 0),
            cor_N_D_pearson   = safe_cor(N, D, "pearson"),
            cor_N_D_spearman  = safe_cor(N, D, "spearman"),
            cor_ds_dc_pearson = safe_cor(ds, dc, "pearson"),
            cor_ds_dc_spearman= safe_cor(ds, dc, "spearman"),
            cor_s_c_pearson   = safe_cor(s_k, c_k, "pearson"),
            cor_s_c_spearman  = safe_cor(s_k, c_k, "spearman"),
            cor_slo_chi_spearman = safe_cor(s_lo, c_k,  "spearman"),
            cor_slo_clo_spearman = safe_cor(s_lo, c_lo, "spearman"),
            pivot_intercept   = pivot_fit(N, D)[1],
            pivot_theta       = pivot_fit(N, D)[2],
            L_star_median     = pivot_fit(N, D)[3],
            L_star_IQR        = pivot_fit(N, D)[4],
            IQR_N             = IQR(N[is.finite(N)]),
            IQR_D             = IQR(D[is.finite(D)]),
            .groups = "drop") %>%
  arrange(cond_code, boundary)
write.csv(corr_cond, out("correlation_cond.csv"), row.names = FALSE)

# =============================================================================
# 4. 요약 텍스트
# =============================================================================
sink(out("reversal_summary.txt"))
tryCatch({
  cat("NRM 문항 1: 채점함수 전치 시 절편모수와 경계모수\n")
  cat("D_k = a1 (s_k - s_{k-1}),  N_k = c_{k-1} - c_k,  b_k = N_k / D_k\n")
  cat(sprintf("반복: 수렴하고 b1~b4 가 계산된 반복 %d개 (제외 없음)\n\n", nrow(rep_df)))

  cat("1) 전체 조건 합계: 채점함수 상태별\n")
  print(as.data.frame(long %>% group_by(boundary, sf_state) %>% state_summary()),
        digits = 3, row.names = FALSE)

  for (bb in unique(long$boundary)) {
    lb <- long %>% filter(boundary == bb)
    cat(sprintf("\n2) %s: 전치 비율 (행 s4, 열 경계모수 평균 IV3)\n", bb))
    print(round(tapply(lb$D <= 0, list(s4 = lb$iv1_sf4, mean = lb$iv3_b_mean), mean), 3))
    rv <- lb %>% filter(D <= 0)
    if (nrow(rv) > 0) {
      cat(sprintf("\n   %s: 전치 반복 중 절편 간격 N < 0 비율\n", bb))
      print(round(tapply(rv$N < 0, list(s4 = rv$iv1_sf4, mean = rv$iv3_b_mean), mean), 3))
      cat(sprintf("\n   %s: 전치 반복 중 b_hat < 0 비율\n", bb))
      print(round(tapply(rv$b_hat < 0, list(s4 = rv$iv1_sf4, mean = rv$iv3_b_mean), mean), 3))
      cat(sprintf("\n   %s: 전치 반복 중 b_hat < 참값 비율\n", bb))
      print(round(tapply(rv$b_hat < rv$b_true, list(s4 = rv$iv1_sf4, mean = rv$iv3_b_mean), mean), 3))
    }
    cb <- corr_cond %>% filter(boundary == bb)
    cat(sprintf("\n   %s: 조건 내 cor(N, D) Spearman 의 조건 평균 (행 s4, 열 IV3)\n", bb))
    print(round(tapply(cb$cor_N_D_spearman, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean, na.rm = TRUE), 2))
    if (bb == "b4") {
      cat("\n   b4: 조건 내 cor(s3, c4) Spearman 의 조건 평균\n")
      print(round(tapply(cb$cor_slo_chi_spearman, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean, na.rm = TRUE), 2))
      cat("\n   b4: 조건 내 cor(s3, c3) Spearman 의 조건 평균\n")
      print(round(tapply(cb$cor_slo_clo_spearman, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean, na.rm = TRUE), 2))
    }
    cat(sprintf("\n   %s: 조건 내 cor(채점함수 간격, 절편 간격) Spearman 의 조건 평균\n", bb))
    print(round(tapply(cb$cor_ds_dc_spearman, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean, na.rm = TRUE), 2))
    cat(sprintf("\n   %s: 피벗 theta* (N ~ D 회귀 기울기) 의 조건 평균, 참고로 참값 b 평균\n", bb))
    print(round(tapply(cb$pivot_theta, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean, na.rm = TRUE), 2))
    print(round(tapply(cb$b_true, list(s4 = cb$iv1_sf4, mean = cb$iv3_b_mean), mean), 2))
  }
}, finally = sink())
message(sprintf("저장 완료: %s, %s, %s", out("reversal_table_cond.csv"),
                out("correlation_cond.csv"), out("reversal_summary.txt")))

# =============================================================================
# 5. 그래프 (전치 그래프 배치: 행 IV1, 열 IV2 좁음 → 넓음, x축 IV3)
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

for (bb in unique(long$boundary)) {
  cd <- add_factors(corr_cond %>% filter(boundary == bb)) %>%
    select(iv1, iv2, iv3, cor_N_D_spearman, cor_ds_dc_spearman) %>%
    pivot_longer(c(cor_N_D_spearman, cor_ds_dc_spearman), names_to = "what", values_to = "r") %>%
    mutate(what = factor(what, levels = c("cor_N_D_spearman", "cor_ds_dc_spearman"),
                         labels = c("cor(N, D): intercept gap vs alpha*scoring gap",
                                    "cor(ds, dc): scoring gap vs intercept gap")))
  p_c <- ggplot(cd, aes(x = iv3, y = r, colour = what, group = what)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_colour_manual(values = c("#4C9BE8", "#D1495B"), name = NULL) +
    scale_y_continuous(limits = c(-1, 1)) +
    labs(title = sprintf("Item 1 NRM %s: within-condition correlation across replications (Spearman)", bb),
         y = "Correlation") + grid_layers
  ggsave(out(sprintf("corr_N_D_%s.png", bb)), p_c, width = 24, height = 18, dpi = 150)

  rv <- add_factors(rev_cond %>% filter(boundary == bb))
  sh <- rv %>% filter(sf_state == "reversed (D <= 0)") %>%
    transmute(iv1, iv2, iv3, `Share of replications reversed (D <= 0)` = share_of_condition,
              `Among reversed: share b_hat < 0` = share_b_negative,
              `Among reversed: share intercept gap N < 0` = share_N_negative) %>%
    pivot_longer(-c(iv1, iv2, iv3), names_to = "what", values_to = "share")
  p_s <- ggplot(sh, aes(x = iv3, y = share, colour = what, group = what)) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_colour_manual(values = c("grey30", "#D1495B", "#4C9BE8"), name = NULL) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = sprintf("Item 1 NRM %s: scoring-function reversal and what happens to the intercept gap and b", bb),
         y = "Share") + grid_layers
  ggsave(out(sprintf("reversed_share_b_neg_%s.png", bb)), p_s, width = 24, height = 18, dpi = 150)
}
message(sprintf("그래프 저장 완료 → %s/", OUT_DIR))
