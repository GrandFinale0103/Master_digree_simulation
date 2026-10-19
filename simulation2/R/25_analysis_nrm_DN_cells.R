# =============================================================================
# 25_analysis_nrm_DN_cells.R  (simulation2, NRM 추정 결과)
# 분모 D(기울기 × 채점함수 간격)와 분자 N(절편 간격)의 추정 상태 조합별 경계모수 편향 (문항 1)
#
# ── 정의 (경계 k, 범주 k 와 k+1 사이) ────────────────────────────────────────
#   D = α (s_k − s_{k−1}),  N = c_{k−1} − c_k,  b = N / D
#   참값: D = sf_k − sf_{k−1} (α = 1),  N = D × b  (둘 다 양수)
#   상태 (D, N 각각):  over (추정 ≥ 참값) / under (0 < 추정 < 참값) / reversed (추정 ≤ 0)
#   → 3 × 3 = 9 칸
#
# ── 상대오차 ────────────────────────────────────────────────────────────────
#   r_D = D̂ / D,  r_N = N̂ / N,  b̂ / b = r_N / r_D
#   D̂ > 0, N̂ > 0 이면 b̂ < b  ⇔  r_N < r_D   (방향은 절댓값이 아니라 "참값 대비 비율" 차이로 정해짐)
#   D̂ ≤ 0, N̂ > 0 이면 b̂ < 0 < b (항상 과소)
#   D̂ ≤ 0, N̂ ≤ 0 이면 b̂ > 0, b̂ < b ⇔ |r_N| < |r_D|
#   N 이 D 를 따라가는 정도: (r_N − 1) = κ (r_D − 1) 회귀 (원점 통과)
#     κ = 1 이면 완전 보정(b 불변), κ < 1 이면 D 가 커질 때 b 가 줄고 D 가 작아질 때 b 가 커짐
#
# ── 출력 (output/nrm/DN_cells/) ───────────────────────────────────────────────
#   DN_cells_cond.csv      — 조건 × 경계 × 9칸: n, 비율, b 중앙값 편향, 참값보다 작은 비율, r_D·r_N 중앙값
#   rD_bins_cond.csv       — 조건 × 경계 × r_D 구간: n, r_N 중앙값, b 중앙값 편향, 참값보다 작은 비율
#   kappa_cond.csv         — 조건 × 경계: κ, R², D>0 중 r_N < r_D 비율, 전체 중앙값 편향, 참값보다 작은 비율
#   DN_cells_summary.txt   — 전체 합계 9칸 표, r_D 구간 표, κ·과소추정 비율의 s4 × IV3 표
#   DN_cells_share_<경계>.png — 조건별 9칸 비율 누적 막대 (색 = 과소추정 칸 / 과대추정 칸 구분)
#
# 반복 집합: 수렴하고 b1~b4 가 계산된 반복 전체 (제외 없음)
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
BOUNDARIES <- 4                                          # 1:4 로 바꾸면 b1~b4 모두
RD_BREAKS  <- c(-Inf, -1, 0, 0.5, 1, 1.5, 2, 3, Inf)     # r_D = D̂/D 구간
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
OUT_DIR <- "output/nrm/DN_cells"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)
STATE_LEVELS <- c("over", "under", "reversed")

# =============================================================================
# 1. 데이터
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

state <- function(est, true) factor(ifelse(est <= 0, "reversed", ifelse(est >= true, "over", "under")),
                                    levels = STATE_LEVELS)
long <- bind_rows(lapply(BOUNDARIES, function(k) {
  sf     <- cbind(0, 1, 2, rep_df$iv1_sf4, 4)
  D_true <- sf[, k + 1] - sf[, k]
  b_true <- rep_df$iv3_b_mean + rep_df$iv2_b_interval * B_OFFSET[k]
  D_hat  <- rep_df$a1 * (rep_df[[paste0("ak", k)]] - rep_df[[paste0("ak", k - 1)]])
  N_hat  <- rep_df[[paste0("d", k - 1)]] - rep_df[[paste0("d", k)]]
  data.frame(cond_code = rep_df$cond_code, iv1_sf4 = rep_df$iv1_sf4,
             iv2_b_interval = rep_df$iv2_b_interval, iv3_b_mean = rep_df$iv3_b_mean,
             boundary = paste0("b", k), b_true = b_true,
             D_true = D_true, N_true = D_true * b_true, D_hat = D_hat, N_hat = N_hat,
             b_hat = N_hat / D_hat, stringsAsFactors = FALSE)
})) %>%
  filter(is.finite(b_hat), N_true > 0) %>%            # N 참값이 양수인 경계만 (b > 0)
  mutate(err = b_hat - b_true, rD = D_hat / D_true, rN = N_hat / N_true,
         D_state = state(D_hat, D_true), N_state = state(N_hat, N_true),
         cell = paste0("D ", D_state, " / N ", N_state),
         below = b_hat < b_true)
dropped <- setdiff(paste0("b", BOUNDARIES), unique(long$boundary))
if (length(dropped) > 0) message("참값 b ≤ 0 인 반복만 있어 제외된 경계: ", paste(dropped, collapse = ", "))

cell_stats <- function(df) df %>% summarise(n = n(), median_bias = median(err), share_below_true = mean(below),
                                            median_rD = median(rD), median_rN = median(rN), .groups = "drop")

# =============================================================================
# 2. 표
# =============================================================================
cells_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, D_state, N_state) %>%
  cell_stats() %>%
  group_by(cond_code, boundary) %>% mutate(share = n / sum(n)) %>% ungroup()
write.csv(cells_cond, out("DN_cells_cond.csv"), row.names = FALSE)

long$rD_bin <- cut(long$rD, RD_BREAKS, right = FALSE)
bins_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, rD_bin) %>%
  summarise(n = n(), median_rN = median(rN), median_b_hat = median(b_hat), median_bias = median(err),
            share_below_true = mean(below), .groups = "drop") %>%
  group_by(cond_code, boundary) %>% mutate(share = n / sum(n)) %>% ungroup()
write.csv(bins_cond, out("rD_bins_cond.csv"), row.names = FALSE)

kappa_fit <- function(rN, rD) {
  ok <- is.finite(rN) & is.finite(rD)
  if (sum(ok) < 10) return(c(NA, NA))
  f <- lm(I(rN[ok] - 1) ~ 0 + I(rD[ok] - 1))
  c(unname(coef(f)), summary(f)$r.squared)
}
kappa_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, b_true, D_true) %>%
  summarise(n = n(),
            kappa = kappa_fit(rN, rD)[1], kappa_R2 = kappa_fit(rN, rD)[2],
            share_rN_lt_rD_Dpos = mean((rN < rD)[D_hat > 0]),
            share_D_over = mean(D_state == "over"), share_D_under = mean(D_state == "under"),
            share_D_reversed = mean(D_state == "reversed"),
            median_bias = median(err), share_below_true = mean(below), .groups = "drop")
write.csv(kappa_cond, out("kappa_cond.csv"), row.names = FALSE)

sink(out("DN_cells_summary.txt"))
tryCatch({
  cat("NRM 문항 1: 분모 D 와 분자 N 의 추정 상태 조합별 경계모수 편향\n")
  cat("D = a1 (s_k - s_{k-1}), N = c_{k-1} - c_k, b = N / D;  b_hat/b = r_N / r_D\n")
  cat(sprintf("반복: %d개 (제외 없음)\n", nrow(rep_df)))
  for (bb in unique(long$boundary)) {
    lb <- long %>% filter(boundary == bb)
    cat(sprintf("\n==== %s ====\n1) 9칸 (전체 조건 합계)\n", bb))
    print(as.data.frame(lb %>% group_by(D_state, N_state) %>% cell_stats() %>%
                          mutate(share = n / sum(n))), digits = 3, row.names = FALSE)
    cat(sprintf("\n   전체: 중앙값 편향 %.3f, 참값보다 작은 비율 %.3f\n", median(lb$err), mean(lb$below)))
    cat("\n2) r_D = D_hat/D 구간별 (전체 조건 합계)\n")
    print(as.data.frame(lb %>% group_by(rD_bin) %>%
                          summarise(n = n(), share = n() / nrow(lb), median_rN = median(rN),
                                    median_bias = median(err), share_below_true = mean(below), .groups = "drop")),
          digits = 3, row.names = FALSE)
    kb <- kappa_cond %>% filter(boundary == bb)
    cat("\n3) kappa: (r_N - 1) = kappa (r_D - 1) 의 조건 평균 (행 s4, 열 IV3)  — 1 이면 완전 보정\n")
    print(round(tapply(kb$kappa, list(s4 = kb$iv1_sf4, mean = kb$iv3_b_mean), mean, na.rm = TRUE), 2))
    cat("\n4) D>0 반복 중 r_N < r_D (b 과소) 비율의 조건 평균\n")
    print(round(tapply(kb$share_rN_lt_rD_Dpos, list(s4 = kb$iv1_sf4, mean = kb$iv3_b_mean), mean, na.rm = TRUE), 3))
    cat("\n5) 전체 반복 중 b_hat < 참값 비율의 조건 평균\n")
    print(round(tapply(kb$share_below_true, list(s4 = kb$iv1_sf4, mean = kb$iv3_b_mean), mean), 3))
    cat("\n6) D 상태 비율의 조건 평균 (over / under / reversed)\n")
    print(round(sapply(c(over = "share_D_over", under = "share_D_under", reversed = "share_D_reversed"),
                       function(v) mean(kb[[v]])), 3))
  }
}, finally = sink())
message(sprintf("저장 완료: %s/ (DN_cells_cond.csv, rD_bins_cond.csv, kappa_cond.csv, DN_cells_summary.txt)", OUT_DIR))

# =============================================================================
# 3. 그래프: 조건별 9칸 비율 (전치 그래프 배치)
# =============================================================================
make_num_factor <- function(x, fmt) {
  vals <- sort(unique(na.omit(x))); factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
CELL_ORDER <- as.vector(outer(paste0("D ", STATE_LEVELS), paste0(" / N ", STATE_LEVELS), paste0))
for (bb in unique(long$boundary)) {
  pd <- cells_cond %>% filter(boundary == bb) %>%
    mutate(cell = factor(paste0("D ", D_state, " / N ", N_state), levels = CELL_ORDER),
           iv1 = make_num_factor(iv1_sf4, "%.2f"), iv2 = make_num_factor(iv2_b_interval, "%.2f"),
           iv3 = make_num_factor(iv3_b_mean, "%.1f"))
  p <- ggplot(pd, aes(x = iv3, y = share, fill = cell)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x), iv2 = function(x) paste0("Int=", x))) +
    scale_fill_manual(values = setNames(c("#08519C", "#6BAED6", "#C6DBEF",
                                          "#FDAE6B", "#E6550D", "#A63603",
                                          "#BCBDDC", "#756BB1", "#3F007D"), CELL_ORDER),
                      name = NULL, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(title = sprintf("Item 1 NRM %s: share of replications by state of D (denominator) and N (numerator)", bb),
         x = "b-parameter Mean (IV3)", y = "Share of replications") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", strip.text = element_text(size = 7),
          axis.text = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(nrow = 3, byrow = FALSE))
  ggsave(out(sprintf("DN_cells_share_%s.png", bb)), p, width = 24, height = 18, dpi = 150)
}
message(sprintf("그래프 저장 완료 → %s/", OUT_DIR))
