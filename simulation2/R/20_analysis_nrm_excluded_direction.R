# =============================================================================
# 20_analysis_nrm_excluded_direction.R  (simulation2, NRM 추정 결과)
# 극단값으로 제외되는 문항 1 경계모수 추정값이 주로 과대추정인지 과소추정인지 확인
#
# ── 1) 극단 추정값의 방향 (경계모수 추정값 단위) ─────────────────────────────
#   극단: |b̂_k − b_k| ≥ ERR_CUT   →  과대(over): b̂_k − b_k ≥ +ERR_CUT,  과소(under): ≤ −ERR_CUT
#   경계모수(b1~b4)별·조건별로 과대/과소 수와 과대 비율
#   부호 검정: H0 과대 비율 = 0.5 (이항검정, 양측), 조건별 검정은 BH(FDR) 보정
#
# ── 2) 반복 제외가 Bias 를 어느 쪽으로 옮기는가 (반복 단위, 16번과 같은 기준) ──
#   16번은 b1~b4 중 하나라도 극단이면 그 반복 전체를 제외하므로,
#   극단이 아니었던 다른 경계모수도 함께 빠짐 → 경계모수마다 제외 효과를 따로 봄
#     shift = Bias(전체) − Bias(제외 후)
#       shift > 0 : 제외된 반복이 평균적으로 과대추정 쪽  → 제외하면 Bias 가 내려감
#       shift < 0 : 제외된 반복이 평균적으로 과소추정 쪽  → 제외하면 Bias 가 올라감
#   제외된 반복의 오차 평균·중앙값, 과대(>0) 비율도 함께 기록
#   (b̂ 가 ±Inf 인 경우: 방향·중앙값·비율에는 부호대로 넣고, 평균·Bias 에서는 뺌)
#   시나리오: excl_extreme / excl_extreme_sfrev (16번과 같음)
#
# ── 3) 방향과 원인 (극단 추정값만) ─────────────────────────────────────────────
#   b_k = (c_{k−1} − c_k) / (α · Δs_k),  Δs_k = s_k − s_{k−1}
#   Δs_k 상태: 전치(Δs_k ≤ 0) / 분모 붕괴(0 < Δs_k < DAK_CUT) / 정상(Δs_k ≥ DAK_CUT)
#   → 방향 × Δs_k 상태 교차표 (분모 부호가 뒤집히면 방향도 뒤집히는지 확인)
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/excluded_direction/
#   direction_summary.txt         — 경계모수별 과대/과소 수·비율·부호 검정, 원인 교차표, 제외 효과 요약
#   extreme_direction_overall.csv — 경계모수별 극단 추정값 방향 (전체 조건 합계)
#   extreme_direction_cond.csv    — 조건 × 경계모수별 극단 추정값 방향 + 부호 검정(BH 보정)
#   extreme_direction_cause.csv   — 경계모수 × 방향 × Δs_k 상태 빈도
#   exclusion_shift_cond.csv      — 시나리오 × 조건 × 경계모수별 Bias(전체/제외 후), shift, 제외 반복의 오차
#   extreme_direction_b<k>.png    — 조건별 과대(위, +)·과소(아래, −) 극단 비율 막대 (경계모수마다 1장)
#   extreme_error_hist.png        — 극단 추정값의 부호 있는 오차 분포 (경계모수별, pseudo-log 축)
#   exclusion_shift_<시나리오>.png  — 제외로 인한 Bias 이동량 shift (경계모수별 색)
#   (조건 그래프 배치: 행 IV1(s4), 열 IV2(간격, 좁음→넓음), x축 IV3(평균) — 전치 그래프와 같음)
# =============================================================================

ERR_CUT <- 3     # 극단 오차: |b̂ − b| ≥ ERR_CUT   (15·16번과 같은 값)
DAK_CUT <- 0.1   # 분모 붕괴: |Δs_k| < DAK_CUT      (15·16번과 같은 값)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
OUT_DIR <- "output/nrm/excluded_direction"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)
B_PARAMS <- paste0("b", 1:4)
DIR_LEVELS  <- c("Over (estimate - true >= cut)", "Under (estimate - true <= -cut)")
DIR_COLOURS <- setNames(c("#D1495B", "#4C9BE8"), DIR_LEVELS)
DAK_LEVELS  <- c("reversed (ds <= 0)", "near-zero (0 < ds < cut)", "normal (ds >= cut)")

# =============================================================================
# 1. 데이터: 반복 × 경계모수 (16번과 같은 반복 집합·제외 기준)
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
             iv3_b_mean = rep_df$iv3_b_mean,
             param = B_PARAMS[k], err = err[, k], dak = dak[, k],
             rep_extreme = rep_df$is_extreme, rep_sfrev = rep_df$is_sfrev,
             stringsAsFactors = FALSE)
})) %>%
  filter(!is.nan(err)) %>%                       # NaN(0/0) 은 방향을 정할 수 없어 제외
  mutate(fin       = is.finite(err),                # ±Inf 는 방향 집계에는 넣고 평균·Bias 에서는 뺌
         extreme   = abs(err) >= ERR_CUT,
         direction = case_when(err >=  ERR_CUT ~ DIR_LEVELS[1],
                               err <= -ERR_CUT ~ DIR_LEVELS[2],
                               TRUE            ~ NA_character_),
         dak_state = case_when(dak <= 0       ~ DAK_LEVELS[1],
                               dak <  DAK_CUT ~ DAK_LEVELS[2],
                               TRUE           ~ DAK_LEVELS[3]))

sign_test <- function(n_over, n_ext) {
  if (n_ext == 0) return(NA_real_)
  binom.test(n_over, n_ext, p = 0.5)$p.value
}

# =============================================================================
# 2. 극단 추정값의 방향
# =============================================================================
dir_overall <- long %>%
  group_by(param) %>%
  summarise(n_estimates  = n(),
            n_extreme    = sum(extreme),
            n_over       = sum(err >=  ERR_CUT),
            n_under      = sum(err <= -ERR_CUT),
            prop_over    = n_over / n_extreme,
            n_infinite   = sum(!fin),
            median_err_extreme = median(err[extreme]),
            .groups = "drop") %>%
  rowwise() %>% mutate(p_sign_test = sign_test(n_over, n_extreme)) %>% ungroup()
write.csv(dir_overall, out("extreme_direction_overall.csv"), row.names = FALSE)

dir_cond <- long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
  summarise(n_converged = n(),
            n_extreme   = sum(extreme),
            n_over      = sum(err >=  ERR_CUT),
            n_under     = sum(err <= -ERR_CUT),
            prop_over_of_extreme = ifelse(n_extreme > 0, n_over / n_extreme, NA_real_),
            rate_over   = n_over / n_converged,
            rate_under  = n_under / n_converged,
            .groups = "drop") %>%
  rowwise() %>% mutate(p_sign_test = sign_test(n_over, n_extreme)) %>% ungroup() %>%
  mutate(p_sign_test_BH = p.adjust(p_sign_test, method = "BH"),
         majority = case_when(n_extreme == 0          ~ "none",
                              p_sign_test_BH >= 0.05  ~ "no clear direction",
                              n_over > n_under        ~ "mostly over",
                              TRUE                    ~ "mostly under"))
write.csv(dir_cond, out("extreme_direction_cond.csv"), row.names = FALSE)

dir_cause <- long %>% filter(extreme) %>%
  count(param, direction, dak_state, name = "n") %>%
  group_by(param, direction) %>% mutate(prop_within_direction = n / sum(n)) %>% ungroup()
write.csv(dir_cause, out("extreme_direction_cause.csv"), row.names = FALSE)

# =============================================================================
# 3. 반복 제외가 Bias 를 옮기는 방향
# =============================================================================
scenarios <- list(
  excl_extreme       = function(d) d$rep_extreme,
  excl_extreme_sfrev = function(d) d$rep_extreme | d$rep_sfrev
)
shift_cond <- bind_rows(lapply(names(scenarios), function(sc) {
  long %>% mutate(excluded = scenarios[[sc]](.)) %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
    summarise(n_converged      = n(),
              n_excluded       = sum(excluded),
              bias_all         = mean(err[fin]),
              bias_kept        = if (any(!excluded & fin)) mean(err[!excluded & fin]) else NA_real_,
              shift            = bias_all - bias_kept,
              mean_err_excluded   = if (any(excluded & fin)) mean(err[excluded & fin]) else NA_real_,
              median_err_excluded = if (any(excluded)) median(err[excluded]) else NA_real_,
              prop_over_excluded  = if (any(excluded)) mean(err[excluded] > 0) else NA_real_,
              .groups = "drop") %>%
    mutate(scenario = sc, .before = 1)
}))
write.csv(shift_cond, out("exclusion_shift_cond.csv"), row.names = FALSE)

shift_overall <- bind_rows(lapply(names(scenarios), function(sc) {
  long %>% mutate(excluded = scenarios[[sc]](.)) %>%
    group_by(param) %>%
    summarise(n_excluded = sum(excluded),
              mean_err_excluded   = mean(err[excluded & fin]),
              median_err_excluded = median(err[excluded]),
              prop_over_excluded  = mean(err[excluded] > 0),
              mean_err_kept       = mean(err[!excluded & fin]),
              .groups = "drop") %>%
    mutate(scenario = sc, .before = 1)
}))

# =============================================================================
# 4. 요약 텍스트
# =============================================================================
sink(out("direction_summary.txt"))
tryCatch({
  cat("NRM 문항 1 경계모수: 극단 추정값의 방향 (과대 vs 과소)\n")
  cat(sprintf("극단 기준 |b_hat - b_true| >= %g, 수렴하고 b1~b4 가 계산된 반복 %d개\n\n",
              ERR_CUT, nrow(rep_df)))
  cat("1) 경계모수별 극단 추정값 방향 (전체 조건 합계, 부호 검정 H0: 과대 비율 = 0.5)\n")
  print(as.data.frame(dir_overall), digits = 4, row.names = FALSE)
  cat("\n   조건별 다수 방향 (BH 보정 p < .05 기준)\n")
  print(with(dir_cond, table(param, majority)))
  cat("\n2) 반복 제외 시 제외된 반복의 오차 (전체 조건 합계)\n")
  cat("   mean_err_excluded > 0 이면 제외된 반복은 평균적으로 과대추정 쪽\n")
  print(as.data.frame(shift_overall), digits = 4, row.names = FALSE)
  cat("\n3) 극단 추정값의 방향 × 그 경계의 ds = s_k - s_{k-1} 상태\n")
  print(as.data.frame(dir_cause), digits = 3, row.names = FALSE)
}, finally = sink())
message(sprintf("저장 완료: %s", out("direction_summary.txt")))

# =============================================================================
# 5. 그래프
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

# ── 조건별 과대(+)·과소(−) 극단 비율 (경계모수마다 1장) ──
dir_plot_df <- add_factors(dir_cond) %>%
  select(iv1, iv2, iv3, param, rate_over, rate_under) %>%
  pivot_longer(c(rate_over, rate_under), names_to = "direction", values_to = "rate") %>%
  mutate(rate = ifelse(direction == "rate_under", -rate, rate),
         direction = factor(direction, levels = c("rate_over", "rate_under"), labels = DIR_LEVELS))
for (pp in B_PARAMS) {
  p <- ggplot(dir_plot_df %>% filter(param == pp), aes(x = iv3, y = rate, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = DIR_COLOURS, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = function(x) scales::percent(abs(x), accuracy = 1)) +
    labs(title = sprintf("Item 1 NRM %s: extreme estimates (|error| >= %g) - over (up) vs under (down)",
                         pp, ERR_CUT),
         y = "Share of converged replications") +
    grid_layers
  ggsave(out(sprintf("extreme_direction_%s.png", pp)), p, width = 24, height = 18, dpi = 150)
}
message(sprintf("저장 완료: %s", out("extreme_direction_b1~b4.png")))

# ── 극단 추정값의 부호 있는 오차 분포 ──
ext <- long %>% filter(extreme, is.finite(err))
if (nrow(ext) > 0) {
  p_hist <- ggplot(ext, aes(x = err, fill = direction)) +
    geom_histogram(bins = 80) +
    geom_vline(xintercept = c(-ERR_CUT, ERR_CUT), linetype = "dashed", colour = "grey40") +
    facet_wrap(~ param, ncol = 2, scales = "free_y") +
    scale_x_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(-1000, -100, -10, -3, 0, 3, 10, 100, 1000)) +
    scale_fill_manual(values = DIR_COLOURS, name = NULL) +
    labs(title = "Signed error of extreme NRM boundary estimates (pseudo-log axis)",
         x = "estimate - true", y = "Count") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(out("extreme_error_hist.png"), p_hist, width = 12, height = 8, dpi = 150)
  message(sprintf("저장 완료: %s", out("extreme_error_hist.png")))
}

# ── 제외로 인한 Bias 이동량 (시나리오마다 1장) ──
param_colours <- c(b1 = "#4C9BE8", b2 = "#59A14F", b3 = "#F2A541", b4 = "#D1495B")
for (sc in names(scenarios)) {
  sd_ <- add_factors(shift_cond %>% filter(scenario == sc))
  p <- ggplot(sd_, aes(x = iv3, y = shift, colour = param, group = param)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.4) +
    scale_colour_manual(values = param_colours, name = NULL) +
    labs(title = sprintf("Bias shift from excluding replications (%s): Bias(all) - Bias(kept); > 0 = excluded were overestimates", sc),
         y = "Bias(all) - Bias(kept)") +
    grid_layers
  ggsave(out(sprintf("exclusion_shift_%s.png", sc)), p, width = 24, height = 18, dpi = 150)
}
message(sprintf("저장 완료: %s", out("exclusion_shift_<scenario>.png")))
