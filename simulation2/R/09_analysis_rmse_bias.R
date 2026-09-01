# =============================================================================
# 09_analysis_rmse_bias.R  (simulation2)
# 문항 1 추정 모수의 편향(Bias)과 RMSE 분석
#
# simulation2 변경:
#   - 조건 코드: 8자리 2자리 고정 형식 (예: "01030201")
#   - IV1: 7 수준, IV2: 10 수준, IV3: 4 수준, IV4: 1 수준 (정규분포만)
#   - 총 280 조건
#
# 출력 파일:
#   output/analysis/09_log.txt
#   output/analysis/rmse_bias_long.csv
#   output/analysis/rmse_bias_summary.csv
#   output/analysis/rmse_bias_plot_b_bias.png
#   output/analysis/rmse_bias_plot_b_rmse.png
#   output/analysis/rmse_bias_plot_a.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── 로그 시스템 ───────────────────────────────────────────────────────────────
LOG_PATH  <- "output/analysis/09_log.txt"
log_lines <- character(0)
lg <- function(...) {
  msg <- paste0(...); cat(msg, "\n"); log_lines <<- c(log_lines, msg)
}
lg_section <- function(title) lg(sprintf("\n[%s]  %s", Sys.time(), title))
flush_log  <- function() writeLines(log_lines, LOG_PATH)

# ── 조건 코드 → 실제 값 매핑 ─────────────────────────────────────────────────
IV1_MAP <- c("01"=3,"02"=3.2,"03"=3.33,"04"=3.4,"05"=3.6,"06"=3.66,"07"=3.8)
IV2_MAP <- c("01"=1.4,"02"=1.25,"03"=1.2,"04"=1,"05"=0.8,
             "06"=0.6,"07"=0.5,"08"=0.4,"09"=0.25,"10"=0.2)
IV3_MAP <- c("01"=0,"02"=0.5,"03"=1,"04"=1.5)
IV4_MAP <- c("01"="normal")
DISCRIM <- 1

EXPECTED_CONDS <- character(0)
for (d1 in seq_along(IV1_MAP)) for (d2 in seq_along(IV2_MAP))
  for (d3 in seq_along(IV3_MAP)) for (d4 in seq_along(IV4_MAP))
    EXPECTED_CONDS <- c(EXPECTED_CONDS, sprintf("%02d%02d%02d%02d", d1, d2, d3, d4))
N_CONDS <- length(EXPECTED_CONDS)  # 280

# ── 파일 검증 함수 ────────────────────────────────────────────────────────────
validate_est_files <- function(files) {
  lg_section("파일 검증 (est_params)")
  lg(sprintf("  발견 파일 수: %d", length(files)))

  parse_cond <- function(p) {
    m <- regmatches(p, regexpr("_cond([0-9]{8})_", p))
    if (length(m) == 0) return(NA_character_)
    gsub("_cond|_", "", m)
  }

  found_conds <- sapply(files, parse_cond)
  reps_tbl    <- sort(table(found_conds[!is.na(found_conds)]))

  missing_conds <- setdiff(EXPECTED_CONDS, names(reps_tbl))
  extra_conds   <- setdiff(names(reps_tbl), EXPECTED_CONDS)

  if (length(missing_conds) > 0) {
    lg(sprintf("  [경고] 누락 조건 %d개", length(missing_conds)))
  } else {
    lg(sprintf("  누락 조건 없음 (%d개 전체 존재)", N_CONDS))
  }
  if (length(extra_conds) > 0)
    lg(sprintf("  [경고] 예상 외 조건 %d개", length(extra_conds)))

  if (length(reps_tbl) > 0) {
    max_reps    <- max(reps_tbl)
    under_conds <- reps_tbl[reps_tbl < max_reps]
    lg(sprintf("  최대 반복 횟수: %d", max_reps))
    if (length(under_conds) > 0)
      lg(sprintf("  [경고] 반복 미달 조건 %d개 (기준: %d회)",
                 length(under_conds), max_reps))
    else
      lg(sprintf("  모든 조건이 동일 반복 횟수(%d회) 충족", max_reps))
  }
  invisible(NULL)
}

# =============================================================================
# 1단계: est_params 파일 로딩 (문항 1만)
# =============================================================================
lg_section("1단계: est_params 파일 로딩")

est_files <- list.files(
  path = "output/estimated_params",
  pattern = "_est_params\\.csv$", full.names = TRUE
)
if (length(est_files) == 0) {
  flush_log()
  stop("추정 결과 파일이 없습니다. 먼저 시뮬레이션을 실행하세요.")
}

validate_est_files(est_files)
lg(sprintf("파일 %d개 로딩 중...", length(est_files)))

# 8자리 조건 코드 파싱
parse_fname <- function(path) {
  fname <- basename(path)
  cond  <- str_match(fname, "_cond([0-9]{8})_")[, 2]
  rep   <- as.integer(str_match(fname, "_rep([0-9]+)_")[, 2])
  list(cond_code = cond, rep_id = rep)
}

read_item1_est <- function(path) {
  meta <- parse_fname(path)
  df   <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(df)) return(NULL)

  if ("success" %in% names(df) && isFALSE(df$success[1])) {
    return(data.frame(
      cond_code = meta$cond_code, rep_id = meta$rep_id,
      converged = FALSE,
      a_est = NA_real_,
      b1_est = NA_real_, b2_est = NA_real_,
      b3_est = NA_real_, b4_est = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  item1 <- df[df$item %in% c("Item1", "item1", 1), ]
  if (nrow(item1) == 0) item1 <- df[1, ]

  data.frame(
    cond_code = meta$cond_code, rep_id = meta$rep_id,
    converged = as.logical(item1$converged[1]),
    a_est  = item1$a[1],
    b1_est = item1$b1[1], b2_est = item1$b2[1],
    b3_est = item1$b3[1], b4_est = item1$b4[1],
    stringsAsFactors = FALSE
  )
}

est_list <- lapply(est_files, read_item1_est)
est_raw  <- do.call(rbind, est_list[!sapply(est_list, is.null)])
lg(sprintf("  총 %d행 로딩 완료", nrow(est_raw)))

# =============================================================================
# 2단계: 조건 코드 분해 (8자리 2자리 고정) 및 참값 계산
# =============================================================================
est_raw <- est_raw %>%
  mutate(
    k1 = substr(cond_code, 1, 2),
    k2 = substr(cond_code, 3, 4),
    k3 = substr(cond_code, 5, 6),
    k4 = substr(cond_code, 7, 8),
    iv1_sf4        = IV1_MAP[k1],
    iv2_b_interval = IV2_MAP[k2],
    iv3_b_mean     = IV3_MAP[k3],
    iv4_theta_dist = IV4_MAP[k4],
    # 문항 1 참값 (결정론적)
    a_true  = DISCRIM,
    b1_true = iv3_b_mean + iv2_b_interval * (-1.5),
    b2_true = iv3_b_mean + iv2_b_interval * (-0.5),
    b3_true = iv3_b_mean + iv2_b_interval *   0.5,
    b4_true = iv3_b_mean + iv2_b_interval *   1.5
  ) %>%
  select(-k1, -k2, -k3, -k4)

# =============================================================================
# 3단계: 오차 계산
# =============================================================================
est_raw <- est_raw %>%
  mutate(
    err_a  = a_est  - a_true,
    err_b1 = b1_est - b1_true,
    err_b2 = b2_est - b2_true,
    err_b3 = b3_est - b3_true,
    err_b4 = b4_est - b4_true
  )

# =============================================================================
# 4단계: Long format 저장
# =============================================================================
long_df <- est_raw %>%
  filter(converged == TRUE) %>%
  select(cond_code, rep_id,
         iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
         a_est, b1_est, b2_est, b3_est, b4_est,
         a_true, b1_true, b2_true, b3_true, b4_true,
         err_a, err_b1, err_b2, err_b3, err_b4) %>%
  arrange(cond_code, rep_id)

write.csv(long_df, "output/analysis/rmse_bias_long.csv", row.names = FALSE)
lg(sprintf("저장 완료: output/analysis/rmse_bias_long.csv (%d행)", nrow(long_df)))

# =============================================================================
# 5단계: 조건 × 모수별 Bias / RMSE 집계
# =============================================================================
group_vars <- c("cond_code", "iv1_sf4", "iv2_b_interval",
                "iv3_b_mean", "iv4_theta_dist")

compute_metrics <- function(df, err_col, param_name) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n      = n(),
      bias   = mean(.data[[err_col]], na.rm = TRUE),
      rmse   = sqrt(mean(.data[[err_col]]^2, na.rm = TRUE)),
      sd_err = sd(.data[[err_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(param = param_name)
}

summary_df <- bind_rows(
  compute_metrics(long_df, "err_a",  "a"),
  compute_metrics(long_df, "err_b1", "b1"),
  compute_metrics(long_df, "err_b2", "b2"),
  compute_metrics(long_df, "err_b3", "b3"),
  compute_metrics(long_df, "err_b4", "b4")
) %>%
  select(cond_code, param, everything()) %>%
  arrange(cond_code, param)

write.csv(summary_df, "output/analysis/rmse_bias_summary.csv", row.names = FALSE)
lg(sprintf("저장 완료: output/analysis/rmse_bias_summary.csv (%d행)", nrow(summary_df)))

# =============================================================================
# 6단계: 시각화
# =============================================================================
make_iv_factor <- function(x, fmt = "%.4f", decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

summary_df <- summary_df %>%
  mutate(
    iv1 = make_iv_factor(iv1_sf4,        fmt = "%.4f", decreasing = FALSE),
    iv2 = make_iv_factor(iv2_b_interval, fmt = "%.4f", decreasing = TRUE),
    iv3 = make_iv_factor(iv3_b_mean,     fmt = "%.2f", decreasing = FALSE)
  )

b_params <- c("b1", "b2", "b3", "b4")
b_labels <- c(b1 = "b₁", b2 = "b₂", b3 = "b₃", b4 = "b₄")

plot_b <- summary_df %>%
  filter(param %in% b_params) %>%
  mutate(param = factor(param, levels = b_params, labels = b_labels))

# Bias
p_bias_b <- ggplot(plot_b, aes(x = iv3, y = bias, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7, colour = "steelblue") +
  geom_point(size = 1.8, colour = "steelblue") +
  facet_grid(
    rows = vars(param),
    cols = vars(iv2, iv1),
    labeller = labeller(
      iv2   = function(x) paste0("Int=", x),
      iv1   = function(x) paste0("s4=", x),
      param = label_value
    )
  ) +
  labs(x = "b-parameter Mean (IV3)", y = "Bias") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave("output/analysis/rmse_bias_plot_b_bias.png", p_bias_b,
       width = 24, height = 10, dpi = 150)
lg("저장 완료: output/analysis/rmse_bias_plot_b_bias.png")

# RMSE
p_rmse_b <- ggplot(plot_b, aes(x = iv3, y = rmse, group = 1)) +
  geom_line(linewidth = 0.7, colour = "firebrick") +
  geom_point(size = 1.8, colour = "firebrick") +
  facet_grid(
    rows = vars(param),
    cols = vars(iv2, iv1),
    labeller = labeller(
      iv2   = function(x) paste0("Int=", x),
      iv1   = function(x) paste0("s4=", x),
      param = label_value
    )
  ) +
  labs(x = "b-parameter Mean (IV3)", y = "RMSE") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 7))

ggsave("output/analysis/rmse_bias_plot_b_rmse.png", p_rmse_b,
       width = 24, height = 10, dpi = 150)
lg("저장 완료: output/analysis/rmse_bias_plot_b_rmse.png")

# a 모수
plot_a <- summary_df %>% filter(param == "a")

p_a <- ggplot(plot_a, aes(x = iv3, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(aes(y = bias),  linewidth = 0.7, colour = "steelblue") +
  geom_point(aes(y = bias), size = 1.8, colour = "steelblue", shape = 16) +
  geom_line(aes(y = rmse),  linewidth = 0.7, colour = "firebrick", linetype = "dotted") +
  geom_point(aes(y = rmse), size = 1.8, colour = "firebrick", shape = 17) +
  facet_grid(
    rows = vars(iv2),
    cols = vars(iv1),
    labeller = labeller(
      iv1 = function(x) paste0("s4=", x),
      iv2 = function(x) paste0("Int=", x)
    )
  ) +
  labs(x = "b-parameter Mean (IV3)", y = "Value",
       caption = "Solid blue = Bias   Dotted red = RMSE") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8))

ggsave("output/analysis/rmse_bias_plot_a.png", p_a,
       width = 14, height = 10, dpi = 150)
lg("저장 완료: output/analysis/rmse_bias_plot_a.png")

# ── 최종 요약 ─────────────────────────────────────────────────────────────────
lg_section("최종 요약")
lg(sprintf("  분석 반복 수 (수렴 성공): %d", nrow(long_df)))
lg(sprintf("  조건 수                 : %d", n_distinct(summary_df$cond_code)))
for (p in c("a", "b1", "b2", "b3", "b4")) {
  sub <- summary_df[summary_df$param == p, ]
  lg(sprintf("  [%s]  Bias: %+.4f ~ %+.4f  |  RMSE: %.4f ~ %.4f",
             p,
             min(sub$bias, na.rm = TRUE), max(sub$bias, na.rm = TRUE),
             min(sub$rmse, na.rm = TRUE), max(sub$rmse, na.rm = TRUE)))
}
flush_log()
cat(sprintf("\n로그 저장 완료: %s\n", LOG_PATH))
