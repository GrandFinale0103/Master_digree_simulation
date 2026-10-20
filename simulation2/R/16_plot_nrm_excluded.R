# =============================================================================
# 16_plot_nrm_excluded.R  (simulation2, NRM 추정 결과)
# 극단 오차 / 채점함수 이상(전치·분모 붕괴) / 경계모수 전치 반복을 제외한 뒤
# 전치 비율, Bias, RMSE 재계산·시각화
#
# 제외는 반복 단위 (전치 비율과 Bias/RMSE가 같은 반복 집합에서 계산되도록)
#   극단 오차      : b1~b4 중 하나라도 |b̂_k − b_k| ≥ ERR_CUT
#   채점함수 이상  : Δak_k = ak_k − ak_{k−1} 이 ≤ 0 (전치) 또는 |Δak_k| < DAK_CUT (분모 붕괴)
#                    인 k 가 하나라도 있음 (k = 1..4)
#   경계모수 전치  : 인접 경계모수 쌍 중 하나라도 b̂_k ≥ b̂_{k+1}
#
# 시나리오
#   excl_extreme        — 1) 극단 오차 반복 제외
#   excl_extreme_sfrev  — 2) 극단 오차 반복 + 채점함수 이상(전치·분모 붕괴) 반복 제외
#   excl_extreme_trans  — 3) 극단 오차 반복 + 경계모수 전치 반복 제외
#   excl_extreme_sfrev_trans — 4) 극단 오차 + 채점함수 이상 + 경계모수 전치 반복 모두 제외
#   ※ 3)·4) 는 전치된 반복을 빼므로 남은 반복의 전치 비율이 0 → 전치 비율 그래프는 만들지 않음
#     (대신 summary_cond.csv 의 prop_transposed 는 0, 제외 비율 그래프에 전치 몫이 표시됨)
#
# 전치 비율 집계 (12·14·15번 NRM 전치 그래프와 같은 기준)
#   전치 = 인접 경계모수 쌍 중 하나라도 b_k ≥ b_{k+1}
#   전치 비율 = (제외되지 않은 반복 중 전치된 반복 수) / (조건의 전체 수렴 반복 수)
#   → 15번 전치 분해 그래프(all_b/transposition_decomposition.png)와 정확히 맞음
#      excl_extreme       = "채점함수 이상만" + "둘 다 아님"
#      excl_extreme_sfrev = "둘 다 아님"
#   Bias/RMSE 는 제외 후 남은 반복으로 계산
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력)
# 출력: output/nrm/excluded/<시나리오>/
#   summary_cond.csv           — 조건별 수렴·남은·제외 반복 수, 전치 비율(분모: 수렴 반복), b1~b4 Bias/RMSE
#   transposition_plot.png     — 전치 비율 (행 IV1, 열 IV2 좁음→넓음, x축 IV3)
#   rmse_bias_plot_b_bias.png  — b1~b4 Bias (12번 NRM 그래프와 같은 배치)
#   rmse_bias_plot_b_rmse.png  — b1~b4 RMSE
#   exclusion_rate_plot.png    — 조건별 제외 비율 (제외 이유 조합별 누적 막대, 전치 그래프 배치)
# 시나리오 비교: output/nrm/excluded/compare/
#   rmse_bias_compare.csv      — 조건 × 경계모수 × 반복 집합(제외 없음 all + 시나리오 4개)의 n, Bias, RMSE
#   rmse_compare_b{k}.png      — 경계모수별 RMSE 를 반복 집합끼리 비교 (행 IV1, 열 IV2 좁음→넓음, x축 IV3)
#   bias_compare_b{k}.png      — 같은 배치의 Bias 비교
# =============================================================================

ERR_CUT <- 3     # 극단 오차: |오차| ≥ ERR_CUT
DAK_CUT <- 0.1   # 분모 붕괴: |Δak| < DAK_CUT

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)

make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4)) %>%
  mutate(
    err_b1 = b1 - (iv3_b_mean + iv2_b_interval * -1.5),
    err_b2 = b2 - (iv3_b_mean + iv2_b_interval * -0.5),
    err_b3 = b3 - (iv3_b_mean + iv2_b_interval *  0.5),
    err_b4 = b4 - (iv3_b_mean + iv2_b_interval *  1.5),
    is_extreme = pmax(abs(err_b1), abs(err_b2), abs(err_b3), abs(err_b4)) >= ERR_CUT,
    is_sfrev   = {
      dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
      rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
    },
    transposed = b1 >= b2 | b2 >= b3 | b3 >= b4
  )

n_all <- rep_df %>% count(cond_code, name = "n_converged")

# flags: 제외 이유별 논리 벡터 목록 (예: list("Extreme error" = is_extreme, ...)) — 하나라도 TRUE 면 제외
run_scenario <- function(flags, label, desc) {
  keep    <- !Reduce(`|`, flags)
  trans_excluded <- "Transposed b" %in% names(flags)
  out_dir <- file.path("output/nrm/excluded", label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)
  kept <- rep_df[keep, ]

  # 분모 = 조건의 전체 수렴 반복 (제외된 반복도 분모에 남김)
  trans_cond <- rep_df %>%
    mutate(kept_rep = keep) %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    summarise(n_converged  = n(),
              n_kept       = sum(kept_rep),
              n_transposed = sum(transposed & kept_rep), .groups = "drop") %>%
    mutate(n_excluded      = n_converged - n_kept,
           prop_excluded   = n_excluded / n_converged,
           prop_transposed = n_transposed / n_converged)

  err_cond <- kept %>%
    select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, err_b1:err_b4) %>%
    pivot_longer(err_b1:err_b4, names_to = "param", names_prefix = "err_", values_to = "err") %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
    summarise(bias = mean(err), rmse = sqrt(mean(err^2)), .groups = "drop")

  summary_cond <- trans_cond %>%
    left_join(err_cond %>%
                pivot_wider(id_cols = cond_code, names_from = param,
                            values_from = c(bias, rmse)),
              by = "cond_code")
  write.csv(summary_cond, out("summary_cond.csv"), row.names = FALSE)

  # ── 전치 비율 (기존 NRM 전치 그래프 배치) — 전치 반복을 제외한 시나리오는 항상 0 이라 생략 ──
  if (!trans_excluded) {
  p_t <- trans_cond %>%
    mutate(iv1 = make_num_factor(iv1_sf4, "%.2f"),
           iv2 = make_num_factor(iv2_b_interval, "%.2f", decreasing = FALSE),
           iv3 = make_num_factor(iv3_b_mean, "%.1f")) %>%
    ggplot(aes(x = iv3, y = prop_transposed)) +
    geom_col(fill = "steelblue") +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = desc, x = "b-parameter Mean (IV3)",
         y = "Reversal Rate (NRM, any adjacent pair; denominator = all converged reps)") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out("transposition_plot.png"), p_t, width = 24, height = 18, dpi = 150)
  }

  # ── Bias / RMSE (12번 NRM 그래프와 같은 배치) ──
  b_params <- c("b1", "b2", "b3", "b4")
  plot_b <- err_cond %>%
    mutate(param = factor(param, levels = b_params,
                          labels = c("b₁", "b₂", "b₃", "b₄")),
           iv1 = make_num_factor(iv1_sf4, "%.2f"),
           iv2 = make_num_factor(iv2_b_interval, "%.2f", decreasing = TRUE),
           iv3 = make_num_factor(iv3_b_mean, "%.1f"))
  common_layers <- list(
    facet_grid(rows = vars(param), cols = vars(iv2, iv1),
               labeller = labeller(iv2 = function(x) paste0("Int=", x),
                                   iv1 = function(x) paste0("s4=", x),
                                   param = label_value)),
    labs(title = desc, x = "b-parameter Mean (IV3)"),
    theme_bw(base_size = 8),
    theme(strip.text = element_text(size = 6), axis.text = element_text(size = 6),
          axis.title = element_text(size = 8),
          axis.text.x = element_text(angle = 45, hjust = 1))
  )
  p_bias <- ggplot(plot_b, aes(x = iv3, y = bias, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_line(linewidth = 0.6, colour = "steelblue") +
    geom_point(size = 1.2, colour = "steelblue") +
    labs(y = "Bias (NRM)") + common_layers
  ggsave(out("rmse_bias_plot_b_bias.png"), p_bias,
         width = 42, height = 12, dpi = 150, limitsize = FALSE)
  p_rmse <- ggplot(plot_b, aes(x = iv3, y = rmse, group = 1)) +
    geom_line(linewidth = 0.6, colour = "firebrick") +
    geom_point(size = 1.2, colour = "firebrick") +
    labs(y = "RMSE (NRM)") + common_layers
  ggsave(out("rmse_bias_plot_b_rmse.png"), p_rmse,
         width = 42, height = 12, dpi = 150, limitsize = FALSE)

  # ── 제외 비율 (제외 이유 조합별 누적 막대, 전치 그래프 배치) ──
  flag_mat <- do.call(cbind, flags)
  reason_of <- apply(flag_mat, 1, function(r) {
    hit <- names(flags)[r]
    if (length(hit) == 0) NA_character_ else paste(hit, collapse = " + ")
  })
  # 조합 순서: 이유 수가 많은 조합 → 적은 조합 (가능한 조합 전부, 범례가 시나리오마다 일관되게)
  combos <- unlist(lapply(rev(seq_along(flags)), function(m)
    apply(combn(names(flags), m), 2, paste, collapse = " + ")))
  EX_LEVELS  <- combos
  ex_colours <- setNames(c("#6A3D9A", "#D1495B", "#F2A541", "#2a9d8f",
                           "#4C9BE8", "#8c564b", "#7f7f7f")[seq_along(EX_LEVELS)], EX_LEVELS)
  ex_cond <- rep_df %>%
    mutate(excluded = !keep, reason = reason_of) %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    mutate(n_converged = n()) %>%
    filter(excluded) %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, n_converged, reason) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prop = n / n_converged,
           reason = factor(reason, levels = EX_LEVELS))
  ex_plot_df <- n_all %>%
    left_join(rep_df %>% distinct(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean), by = "cond_code") %>%
    left_join(ex_cond %>% select(cond_code, reason, prop), by = "cond_code") %>%
    mutate(prop = coalesce(prop, 0),
           iv1 = make_num_factor(iv1_sf4, "%.2f"),
           iv2 = make_num_factor(iv2_b_interval, "%.2f", decreasing = FALSE),
           iv3 = make_num_factor(iv3_b_mean, "%.1f"))
  p_ex <- ggplot(ex_plot_df %>% filter(!is.na(reason)),
                 aes(x = iv3, y = prop, fill = reason)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    geom_blank(data = ex_plot_df, aes(x = iv3), inherit.aes = FALSE) +
    facet_grid(rows = vars(iv1), cols = vars(iv2), drop = FALSE,
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_fill_manual(values = ex_colours, name = "Excluded because", drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(title = desc, x = "b-parameter Mean (IV3)", y = "Excluded replications (%)") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out("exclusion_rate_plot.png"), p_ex, width = 24, height = 18, dpi = 150)

  message(sprintf("[%s] 남은 반복 %d / %d (제외 %.1f%%) → %s/",
                  label, sum(keep), length(keep), 100 * mean(!keep), out_dir))
  invisible(err_cond %>% left_join(kept %>% count(cond_code, name = "n"), by = "cond_code") %>%
              mutate(set = label))
}

F_EXT   <- list("Extreme error" = rep_df$is_extreme)
F_SF    <- list("SF reversal/near-zero dak" = rep_df$is_sfrev)
F_TRANS <- list("Transposed b" = rep_df$transposed)

res <- list(
  run_scenario(F_EXT, "excl_extreme",
               sprintf("Excluding extreme errors (|error| >= %g)", ERR_CUT)),
  run_scenario(c(F_EXT, F_SF), "excl_extreme_sfrev",
               sprintf("Excluding extreme errors (|error| >= %g) and scoring-function reversals / near-zero dak (|dak| < %g)",
                       ERR_CUT, DAK_CUT)),
  run_scenario(c(F_EXT, F_TRANS), "excl_extreme_trans",
               sprintf("Excluding extreme errors (|error| >= %g) and transposed boundaries (any adjacent pair)", ERR_CUT)),
  run_scenario(c(F_EXT, F_SF, F_TRANS), "excl_extreme_sfrev_trans",
               sprintf("Excluding extreme errors (|error| >= %g), scoring-function reversals / near-zero dak (|dak| < %g) and transposed boundaries",
                       ERR_CUT, DAK_CUT))
)

# =============================================================================
# 시나리오 비교: 제외 없음(all) + 시나리오 4개의 Bias / RMSE
# =============================================================================
CMP_DIR <- "output/nrm/excluded/compare"
dir.create(CMP_DIR, recursive = TRUE, showWarnings = FALSE)
SET_LEVELS <- c("all", "excl_extreme", "excl_extreme_sfrev", "excl_extreme_trans", "excl_extreme_sfrev_trans")
SET_LABELS <- c("all (no exclusion)", "extreme", "extreme + SF", "extreme + transposed",
                "extreme + SF + transposed")

all_set <- rep_df %>%
  select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, err_b1:err_b4) %>%
  pivot_longer(err_b1:err_b4, names_to = "param", names_prefix = "err_", values_to = "err") %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
  summarise(n = n(), bias = mean(err), rmse = sqrt(mean(err^2)), .groups = "drop") %>%
  mutate(set = "all")

cmp <- bind_rows(all_set, res) %>%
  left_join(n_all, by = "cond_code") %>%
  mutate(prop_kept = n / n_converged,
         set = factor(set, levels = SET_LEVELS)) %>%
  select(set, cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param,
         n, n_converged, prop_kept, bias, rmse) %>%
  arrange(param, cond_code, set)
write.csv(cmp, file.path(CMP_DIR, "rmse_bias_compare.csv"), row.names = FALSE)

set_colours <- setNames(c("grey35", "#D1495B", "#F2A541", "#2a9d8f", "#4C9BE8"), SET_LABELS)
cmp_plot <- cmp %>%
  mutate(set = factor(SET_LABELS[as.integer(set)], levels = SET_LABELS),
         iv1 = make_num_factor(iv1_sf4, "%.2f"),
         iv2 = make_num_factor(iv2_b_interval, "%.2f", decreasing = FALSE))
for (k in c("b1", "b2", "b3", "b4")) {
  d <- cmp_plot %>% filter(param == k)
  for (what in c("rmse", "bias")) {
    p <- ggplot(d, aes(x = iv3_b_mean, y = .data[[what]], colour = set)) +
      { if (what == "bias") geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) } +
      geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
      facet_grid(rows = vars(iv1), cols = vars(iv2),
                 labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                     iv2 = function(x) paste0("Int=", x))) +
      scale_colour_manual(values = set_colours, name = "Replications") +
      scale_x_continuous(breaks = sort(unique(cmp$iv3_b_mean))) +
      labs(title = sprintf("NRM item 1 %s: %s by exclusion rule", k, toupper(what)),
           x = "b-parameter Mean (IV3)", y = sprintf("%s (NRM)", toupper(what))) +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom",
            strip.text = element_text(size = 7), axis.text = element_text(size = 7),
            axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(CMP_DIR, sprintf("%s_compare_%s.png", what, k)), p,
           width = 24, height = 18, dpi = 150)
  }
}
message(sprintf("시나리오 비교 → %s/", CMP_DIR))
