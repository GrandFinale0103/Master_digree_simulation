# =============================================================================
# 24_analysis_nrm_gap_groups.R  (simulation2, NRM 추정 결과)
# 채점함수 간격 D 의 추정 상태(과대 / 과소 / 전치)별로 반복을 나눠
# 각 그룹의 비율과 경계모수 편향을 정리 (문항 1)
#
# ── 그룹 정의 (경계 k, 범주 k 와 k+1 사이, 1부터) ────────────────────────────
#   D̂ = 추정 채점함수 간격, D = 참 채점함수 간격
#     D_DEF = "ds"       : D = s_k − s_{k−1}            (참값 sf_k − sf_{k−1}, 기본)
#     D_DEF = "alpha_ds" : D = α (s_k − s_{k−1})        (경계모수 식의 분모, 참값은 α = 1 이라 같음)
#   over      : D̂ ≥ D           — 간격 과대   → b̂ 는 받침점 θ* 쪽으로 줄어 과소추정 경향
#   under     : 0 < D̂ < D       — 간격 과소   → b̂ 가 커져 과대추정 경향 (D̂ → 0 이면 매우 큼)
#   reversed  : D̂ ≤ 0           — 채점함수 전치 → b̂ 는 θ* 아래로 떨어져 과소추정, 음수 가능
#   NEAR_TOL > 0 이면 |D̂ − D| / D < NEAR_TOL 인 반복을 near 그룹으로 따로 뺌 (기본 0 = 사용 안 함)
#
# ── 그룹별 지표 (조건 × 경계 × 그룹) ─────────────────────────────────────────
#   n, share (조건 안 비율)
#   D̂ 중앙값, 절편 간격 N̂ = c_{k−1} − c_k 중앙값
#   b̂ 중앙값, 중앙값 편향, 평균 편향(유한값만), b̂ < 참값 비율, b̂ < 0 비율
#   contrib_mean_bias = share × (그룹 평균 편향)  — 그룹들의 합 = 조건의 평균 편향 (유한값 기준)
#   ※ 중앙값 편향은 그룹별로 더해지지 않으므로 그룹 비율과 그룹별 중앙값을 함께 해석
#
# 반복 집합: 수렴하고 b1~b4 가 계산된 반복 전체 (제외 없음)
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/gap_groups/
#   gap_group_cond.csv      — 조건 × 경계 × 그룹별 지표
#   gap_group_overall.csv   — 경계 × 그룹별 전체 조건 합계 + s4 × IV3 별 합계
#   gap_group_summary.txt   — 위 표의 요약 (s4 × IV3 표: 그룹 비율, 그룹별 중앙값 편향)
#   gap_group_share_<경계>.png — 조건별 그룹 비율 누적 막대 (전치 그래프 배치)
#   gap_group_bias_<경계>.png  — 조건별 그룹별 b 중앙값 편향 (같은 배치)
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
BOUNDARIES <- 4            # 요약할 경계모수 (1:4 로 바꾸면 b1~b4 모두)
D_DEF      <- "ds"         # "ds" = 채점함수 간격, "alpha_ds" = α × 채점함수 간격 (경계모수 식 분모)
NEAR_TOL   <- 0            # > 0 이면 |D̂ − D|/D < NEAR_TOL 를 near 그룹으로 분리 (예: 0.1)
BIAS_YLIM  <- c(-6, 6)     # 그룹별 편향 그래프 y축 (범위 밖 값은 가장자리에 표시)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
OUT_DIR <- "output/nrm/gap_groups"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

GROUP_LEVELS <- if (NEAR_TOL > 0) c("over", "near", "under", "reversed") else c("over", "under", "reversed")
GROUP_LABELS <- c(over = "Over (D_hat >= D)", near = sprintf("Near (|D_hat - D|/D < %g)", NEAR_TOL),
                  under = "Under, not reversed (0 < D_hat < D)", reversed = "Reversed (D_hat <= 0)")
GROUP_COLOURS <- c(over = "#4C9BE8", near = "grey60", under = "#F2A541", reversed = "#D1495B")

# =============================================================================
# 1. 데이터: 반복 × 경계
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

long <- bind_rows(lapply(BOUNDARIES, function(k) {
  sf_true <- cbind(0, 1, 2, rep_df$iv1_sf4, 4)
  ds_true <- sf_true[, k + 1] - sf_true[, k]                       # 참 채점함수 간격 (α = 1)
  ds_hat  <- rep_df[[paste0("ak", k)]] - rep_df[[paste0("ak", k - 1)]]
  D_hat   <- if (D_DEF == "alpha_ds") rep_df$a1 * ds_hat else ds_hat
  data.frame(cond_code = rep_df$cond_code, rep_id = rep_df$rep_id,
             iv1_sf4 = rep_df$iv1_sf4, iv2_b_interval = rep_df$iv2_b_interval,
             iv3_b_mean = rep_df$iv3_b_mean, boundary = paste0("b", k),
             D_true = ds_true, D_hat = D_hat,
             N_hat  = rep_df[[paste0("d", k - 1)]] - rep_df[[paste0("d", k)]],
             b_true = rep_df$iv3_b_mean + rep_df$iv2_b_interval * B_OFFSET[k],
             b_hat  = rep_df[[paste0("b", k)]],
             stringsAsFactors = FALSE)
})) %>%
  filter(!is.nan(b_hat)) %>%
  mutate(err   = b_hat - b_true,
         group = case_when(D_hat <= 0                                         ~ "reversed",
                           NEAR_TOL > 0 & abs(D_hat - D_true) / D_true < NEAR_TOL ~ "near",
                           D_hat >= D_true                                    ~ "over",
                           TRUE                                               ~ "under"),
         group = factor(group, levels = GROUP_LEVELS))

# =============================================================================
# 2. 그룹별 지표
# =============================================================================
group_stats <- function(df) df %>%
  summarise(n                  = n(),
            median_D_hat       = median(D_hat),
            median_N_hat       = median(N_hat),
            median_b_hat       = median(b_hat),
            median_bias        = median(err),
            mean_bias          = mean(err[is.finite(err)]),
            share_below_true   = mean(b_hat < b_true),
            share_b_negative   = mean(b_hat < 0),
            .groups = "drop")

cond_tbl <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, D_true, b_true, group, .drop = FALSE) %>%
  group_stats() %>%
  group_by(cond_code, boundary) %>%
  mutate(share = n / sum(n),
         contrib_mean_bias = ifelse(n > 0, share * mean_bias, 0)) %>%
  ungroup() %>%
  filter(!is.na(cond_code)) %>%
  arrange(cond_code, boundary, group)

# 조건 전체 값 (그룹 합) 도 같이 기록
cond_all <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, D_true, b_true) %>%
  group_stats() %>%
  mutate(group = "ALL", share = 1, contrib_mean_bias = mean_bias)
cond_out <- bind_rows(cond_tbl %>% mutate(group = as.character(group)), cond_all) %>%
  mutate(group = factor(group, levels = c(GROUP_LEVELS, "ALL"))) %>%
  arrange(cond_code, boundary, group)
write.csv(cond_out, out("gap_group_cond.csv"), row.names = FALSE)

overall <- bind_rows(
  long %>% group_by(boundary, group) %>% group_stats() %>%
    group_by(boundary) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
    mutate(iv1_sf4 = NA_real_, iv3_b_mean = NA_real_, level = "all conditions"),
  long %>% group_by(boundary, iv1_sf4, iv3_b_mean, group) %>% group_stats() %>%
    group_by(boundary, iv1_sf4, iv3_b_mean) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
    mutate(level = "s4 x IV3"))
write.csv(overall, out("gap_group_overall.csv"), row.names = FALSE)

# =============================================================================
# 3. 요약 텍스트
# =============================================================================
sink(out("gap_group_summary.txt"))
tryCatch({
  cat("NRM 문항 1: 채점함수 간격 D 의 추정 상태별 경계모수 편향\n")
  cat(sprintf("D 정의: %s,  그룹: %s\n", if (D_DEF == "alpha_ds") "alpha * (s_k - s_{k-1})" else "s_k - s_{k-1}",
              paste(GROUP_LABELS[GROUP_LEVELS], collapse = " / ")))
  cat(sprintf("반복: 수렴하고 b1~b4 가 계산된 반복 %d개 (제외 없음)\n\n", nrow(rep_df)))
  cat("1) 전체 조건 합계\n")
  print(as.data.frame(overall %>% filter(level == "all conditions") %>%
                        select(boundary, group, n, share, median_D_hat, median_N_hat,
                               median_bias, mean_bias, share_below_true, share_b_negative)),
        digits = 3, row.names = FALSE)
  for (bb in unique(long$boundary)) {
    ob <- overall %>% filter(level == "s4 x IV3", boundary == bb)
    for (g in GROUP_LEVELS) {
      og <- ob %>% filter(group == g)
      if (nrow(og) == 0) next
      cat(sprintf("\n2) %s, 그룹 %s: 비율 (행 s4, 열 경계모수 평균 IV3)\n", bb, GROUP_LABELS[[g]]))
      print(round(xtabs(share ~ iv1_sf4 + iv3_b_mean, og), 3))
      cat(sprintf("   %s, 그룹 %s: b 중앙값 편향\n", bb, GROUP_LABELS[[g]]))
      print(round(xtabs(median_bias ~ iv1_sf4 + iv3_b_mean, og), 2))
    }
    cat(sprintf("\n3) %s, 전체(그룹 합) 중앙값 편향 (참고)\n", bb))
    ab <- long %>% filter(boundary == bb) %>% group_by(iv1_sf4, iv3_b_mean) %>%
      summarise(m = median(err), .groups = "drop")
    print(round(xtabs(m ~ iv1_sf4 + iv3_b_mean, ab), 2))
  }
}, finally = sink())
message(sprintf("저장 완료: %s, %s, %s", out("gap_group_cond.csv"),
                out("gap_group_overall.csv"), out("gap_group_summary.txt")))

# =============================================================================
# 4. 그래프 (전치 그래프 배치: 행 IV1, 열 IV2 좁음 → 넓음, x축 IV3)
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
fill_scale   <- scale_fill_manual(values = GROUP_COLOURS[GROUP_LEVELS], labels = GROUP_LABELS[GROUP_LEVELS],
                                  name = "Estimated scoring gap", drop = FALSE)
colour_scale <- scale_colour_manual(values = c(GROUP_COLOURS[GROUP_LEVELS], ALL = "black"),
                                    labels = c(GROUP_LABELS[GROUP_LEVELS], ALL = "All replications"),
                                    name = NULL)

for (bb in unique(long$boundary)) {
  sd_ <- add_factors(cond_tbl %>% filter(boundary == bb))
  p_s <- ggplot(sd_, aes(x = iv3, y = share, fill = group)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    fill_scale +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Item 1 NRM %s: share of replications by estimated scoring gap (D = %s)",
                         bb, D_DEF), y = "Share of replications") +
    grid_layers
  ggsave(out(sprintf("gap_group_share_%s.png", bb)), p_s, width = 24, height = 18, dpi = 150)

  bd <- add_factors(cond_out %>% filter(boundary == bb, n > 0)) %>%
    mutate(group = factor(as.character(group), levels = c(GROUP_LEVELS, "ALL")),
           shown = pmin(pmax(median_bias, BIAS_YLIM[1]), BIAS_YLIM[2]))
  p_b <- ggplot(bd, aes(x = iv3, y = shown, colour = group, group = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
    colour_scale +
    coord_cartesian(ylim = BIAS_YLIM) +
    labs(title = sprintf("Item 1 NRM %s: median bias of b by estimated scoring-gap group (values beyond y-limits drawn at the edge)", bb),
         y = "Median bias of b (estimate - true)") +
    grid_layers
  ggsave(out(sprintf("gap_group_bias_%s.png", bb)), p_b, width = 24, height = 18, dpi = 150)
}
message(sprintf("그래프 저장 완료 → %s/", OUT_DIR))
