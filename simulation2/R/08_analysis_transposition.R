# =============================================================================
# 08_analysis_transposition.R  (simulation2)
# 경계 모수 전치(boundary parameter reversal) 분석
#
# simulation2 변경:
#   - 조건 코드: 8자리 2자리 고정 형식 (예: "01030201")
#   - IV1: 7 수준, IV2: 10 수준, IV3: 4 수준, IV4: 1 수준 (정규분포만)
#   - 총 280 조건
#
# 출력 파일:
#   output/analysis/08_log.txt
#   output/analysis/transposition_long.csv
#   output/analysis/transposition_summary.csv
#   output/analysis/transposition_glm.txt
#   output/analysis/transposition_plot.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── 로그 시스템 ───────────────────────────────────────────────────────────────
LOG_PATH  <- "output/analysis/08_log.txt"
log_lines <- character(0)
lg <- function(...) {
  msg <- paste0(...); cat(msg, "\n"); log_lines <<- c(log_lines, msg)
}
lg_section <- function(title) lg(sprintf("\n[%s]  %s", Sys.time(), title))
flush_log  <- function() writeLines(log_lines, LOG_PATH)

# ── 조건 코드 → 실제 값 매핑 ─────────────────────────────────────────────────
# 키: 2자리 인덱스 문자열
IV1_MAP <- c("01"=3,"02"=3.2,"03"=3.33,"04"=3.4,"05"=3.6,"06"=3.66,"07"=3.8)
IV2_MAP <- c("01"=1.4,"02"=1.25,"03"=1.2,"04"=1,"05"=0.8,
             "06"=0.6,"07"=0.5,"08"=0.4,"09"=0.25,"10"=0.2)
IV3_MAP <- c("01"=0,"02"=0.5,"03"=1,"04"=1.5)
IV4_MAP <- c("01"="normal")

# 기대 조건 코드 전체 목록
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
    for (chunk in split(missing_conds, ceiling(seq_along(missing_conds)/10)))
      lg("    ", paste(chunk, collapse = " "))
  } else {
    lg(sprintf("  누락 조건 없음 (%d개 전체 존재)", N_CONDS))
  }
  if (length(extra_conds) > 0)
    lg(sprintf("  [경고] 예상 외 조건 %d개: %s",
               length(extra_conds), paste(extra_conds, collapse = ", ")))

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
# 1단계: est_params 파일 로딩
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

read_est_file <- function(path) {
  meta <- parse_fname(path)
  df   <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(df)) return(NULL)

  if ("success" %in% names(df) && isFALSE(df$success[1])) {
    return(data.frame(
      cond_code = meta$cond_code, rep_id = meta$rep_id,
      converged = FALSE,
      b1 = NA_real_, b2 = NA_real_, b3 = NA_real_, b4 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  item1 <- df[df$item %in% c("Item1", "item1", 1), ]
  if (nrow(item1) == 0) item1 <- df[1, ]

  data.frame(
    cond_code = meta$cond_code, rep_id = meta$rep_id,
    converged = as.logical(item1$converged[1]),
    b1 = item1$b1[1], b2 = item1$b2[1],
    b3 = item1$b3[1], b4 = item1$b4[1],
    stringsAsFactors = FALSE
  )
}

raw_list <- lapply(est_files, read_est_file)
raw      <- do.call(rbind, raw_list[!sapply(raw_list, is.null)])
lg(sprintf("  총 %d행 로딩 완료", nrow(raw)))

# =============================================================================
# 2단계: 조건 코드 분해 (8자리 2자리 고정)
# =============================================================================
lg_section("2단계: 조건 코드 분해")

raw <- raw %>%
  mutate(
    k1 = substr(cond_code, 1, 2),
    k2 = substr(cond_code, 3, 4),
    k3 = substr(cond_code, 5, 6),
    k4 = substr(cond_code, 7, 8),
    iv1_sf4        = IV1_MAP[k1],
    iv2_b_interval = IV2_MAP[k2],
    iv3_b_mean     = IV3_MAP[k3],
    iv4_theta_dist = IV4_MAP[k4]
  ) %>%
  select(-k1, -k2, -k3, -k4)

# =============================================================================
# 3단계: 경계 모수 전치 탐지
# =============================================================================
detect_transposition <- function(b1, b2, b3, b4) {
  if (any(is.na(c(b1, b2, b3, b4)))) {
    return(list(transposed = NA_integer_, which_transposed = NA_character_))
  }
  pairs    <- list(c(1,2), c(2,3), c(3,4))
  bvals    <- c(b1, b2, b3, b4)
  violated <- sapply(pairs, function(p) bvals[p[1]] >= bvals[p[2]])
  which_str <- paste0(
    sapply(which(violated), function(i) paste0(pairs[[i]], collapse = "")),
    collapse = ""
  )
  list(
    transposed       = as.integer(any(violated)),
    which_transposed = if (any(violated)) which_str else ""
  )
}

trans_results        <- mapply(detect_transposition,
                               raw$b1, raw$b2, raw$b3, raw$b4, SIMPLIFY = FALSE)
raw$transposed       <- sapply(trans_results, `[[`, "transposed")
raw$which_transposed <- sapply(trans_results, `[[`, "which_transposed")

# =============================================================================
# 4단계: Long format CSV 저장
# =============================================================================
long_df <- raw %>%
  select(cond_code, rep_id,
         iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
         converged, b1, b2, b3, b4, transposed, which_transposed) %>%
  arrange(cond_code, rep_id)

write.csv(long_df, "output/analysis/transposition_long.csv", row.names = FALSE)
lg(sprintf("저장 완료: output/analysis/transposition_long.csv (%d행)", nrow(long_df)))

# =============================================================================
# 5단계: 조건별 집계
# =============================================================================
summary_df <- long_df %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist) %>%
  summarise(
    n_total      = n(),
    n_converged  = sum(converged,  na.rm = TRUE),
    n_transposed = sum(transposed, na.rm = TRUE),
    prop_transposed = n_transposed / n_converged,
    n_12 = sum(grepl("12", which_transposed), na.rm = TRUE),
    n_23 = sum(grepl("23", which_transposed), na.rm = TRUE),
    n_34 = sum(grepl("34", which_transposed), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, "output/analysis/transposition_summary.csv", row.names = FALSE)
lg(sprintf("저장 완료: output/analysis/transposition_summary.csv (%d행)", nrow(summary_df)))

# =============================================================================
# 6단계: GLM 분석
# =============================================================================
make_iv_factor <- function(x, fmt = "%.2f", decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

IV_CONFIG <- list(
  iv1 = list(col = "iv1_sf4",        type = "numeric", fmt = "%.4f", decreasing = FALSE),
  iv2 = list(col = "iv2_b_interval", type = "numeric", fmt = "%.4f", decreasing = TRUE),
  iv3 = list(col = "iv3_b_mean",     type = "numeric", fmt = "%.2f", decreasing = FALSE),
  iv4 = list(col = "iv4_theta_dist", type = "numeric", fmt = "%.4f", decreasing = FALSE)
)

glm_data <- long_df %>%
  filter(converged == TRUE, !is.na(transposed)) %>%
  mutate(transposed = as.numeric(transposed))

for (iv_name in names(IV_CONFIG)) {
  cfg <- IV_CONFIG[[iv_name]]
  glm_data[[iv_name]] <- make_iv_factor(
    glm_data[[cfg$col]], fmt = cfg$fmt, decreasing = cfg$decreasing
  )
}

iv_vars <- names(IV_CONFIG)

lg_section("6단계: GLM 분석")
lg(sprintf("GLM 분석 대상: %d 반복 (수렴 성공 기준)", nrow(glm_data)))

if (nrow(glm_data) == 0) {
  flush_log()
  stop("수렴 성공 데이터가 없습니다.")
}

lg(sprintf("  transposed 분포: 0=%d건, 1=%d건, NA=%d건",
           sum(glm_data$transposed == 0, na.rm = TRUE),
           sum(glm_data$transposed == 1, na.rm = TRUE),
           sum(is.na(glm_data$transposed))))

if (isTRUE(var(glm_data$transposed, na.rm = TRUE) == 0)) {
  flush_log()
  stop("transposed 분산이 0입니다 — GLM 적합 불가.")
}

lg("  IV별 실제 관측 수준:")
for (v in iv_vars) {
  lg(sprintf("    %s: %d개 수준", v, nlevels(glm_data[[v]])))
}

# 완전 분리 진단
cell_check <- glm_data %>%
  group_by(across(all_of(iv_vars))) %>%
  summarise(n = n(), prop_trans = mean(transposed, na.rm = TRUE), .groups = "drop")
n_sep_zero <- sum(cell_check$prop_trans == 0,   na.rm = TRUE)
n_sep_one  <- sum(cell_check$prop_trans == 1,   na.rm = TRUE)
lg(sprintf("  완전 분리: 0%% 셀=%d개, 100%% 셀=%d개 (전체 %d셀)",
           n_sep_zero, n_sep_one, nrow(cell_check)))

multi_ivs  <- iv_vars[sapply(iv_vars, function(v) nlevels(glm_data[[v]]) >= 2)]
single_ivs <- setdiff(iv_vars, multi_ivs)

if (length(single_ivs) > 0)
  lg(sprintf("  단일 수준 요인 제외: %s", paste(single_ivs, collapse = ", ")))

if (length(multi_ivs) == 0) {
  flush_log()
  stop("GLM에 포함 가능한 요인이 없습니다.")
}

glm_formula <- as.formula(paste("transposed ~", paste(multi_ivs, collapse = " * ")))
lg(sprintf("  GLM 공식: %s", deparse(glm_formula)))

glm_full <- withCallingHandlers(
  glm(glm_formula, data = glm_data, family = binomial(link = "logit")),
  warning = function(w) {
    lg(sprintf("  [GLM 경고] %s", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
)

anova_tbl <- drop1(glm_full, scope = ~., test = "Chisq")
lrt_col   <- intersect(c("LRT", "Chisq"), colnames(anova_tbl))[1]
pval_col  <- grep("^Pr", colnames(anova_tbl), value = TRUE)[1]
lrt_terms <- rownames(anova_tbl)[rownames(anova_tbl) != "<none>"]
partial_r2 <- setNames(
  anova_tbl[lrt_terms, lrt_col] / glm_full$null.deviance,
  lrt_terms
)

# emmeans 대비
if (!requireNamespace("emmeans", quietly = TRUE))
  install.packages("emmeans", repos = "https://cran.rstudio.com/")

emm_list      <- list()
contrast_list <- list()
for (v in multi_ivs) {
  emm_list[[v]]      <- emmeans::emmeans(glm_full,
                                         as.formula(paste("~", v)),
                                         type = "response")
  contrast_list[[v]] <- emmeans::contrast(emm_list[[v]],
                                          method = "pairwise",
                                          adjust = "bonferroni")
}

# GLM 결과 저장
sink("output/analysis/transposition_glm.txt")
cat("================================================================\n")
cat("GLM: boundary parameter reversal ~ condition (simulation2)\n")
cat("family = binomial(logit)\n")
cat("================================================================\n\n")
cat(sprintf("전체: %d반복  수렴 성공: %d반복  전치 발생: %d반복 (%.1f%%)\n\n",
            nrow(long_df),
            sum(long_df$converged, na.rm=TRUE),
            sum(long_df$transposed, na.rm=TRUE),
            100 * mean(long_df$transposed, na.rm=TRUE)))
cat(sprintf("Null deviance: %.2f  Residual deviance: %.2f  McFadden pseudo-R²: %.4f\n\n",
            glm_full$null.deviance, deviance(glm_full),
            1 - deviance(glm_full) / glm_full$null.deviance))
cat("── Type III LRT (drop1) + 편부분 pseudo-R² ──\n\n")
anova_out <- as.data.frame(anova_tbl[lrt_terms, , drop = FALSE])
anova_out$partial_pseudo_R2 <- partial_r2[lrt_terms]
print(anova_out, digits = 4)
cat("\n── Odds Ratio (Wald) ──\n\n")
or_ci <- exp(cbind(OR = coef(glm_full), confint.default(glm_full)))
print(round(or_ci, 4))
cat("\n── emmeans 대비 (Bonferroni) ──\n\n")
iv_labels <- c(iv1 = "IV1 (sf4)", iv2 = "IV2 (b_interval)",
               iv3 = "IV3 (b_mean)", iv4 = "IV4 (dist)")
for (v in multi_ivs) {
  cat(sprintf("\n%s 주변 평균:\n", iv_labels[v]))
  print(emm_list[[v]])
  cat(sprintf("\n%s 쌍별 대비:\n", iv_labels[v]))
  print(contrast_list[[v]])
}
sink()
lg("저장 완료: output/analysis/transposition_glm.txt")

# =============================================================================
# 7단계: 시각화
# =============================================================================
p <- ggplot(summary_df,
            aes(x    = factor(iv3_b_mean),
                y    = prop_transposed,
                fill = factor(iv1_sf4))) +
  geom_col(position = "dodge") +
  facet_grid(
    rows = vars(iv1_sf4),
    cols = vars(iv2_b_interval),
    labeller = labeller(
      iv1_sf4        = function(x) paste0("s4 = ", x),
      iv2_b_interval = function(x) paste0("Interval = ", x)
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "b-parameter Mean (IV3)", y = "Reversal Rate", fill = "s4 (IV1)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/analysis/transposition_plot.png", p,
       width = 12, height = 14, dpi = 150)
lg("저장 완료: output/analysis/transposition_plot.png")

# ── 최종 요약 ─────────────────────────────────────────────────────────────────
lg_section("최종 요약")
lg(sprintf("  전체 반복 수        : %d", nrow(long_df)))
lg(sprintf("  수렴 성공           : %d (%.1f%%)",
           sum(long_df$converged, na.rm=TRUE),
           100 * mean(long_df$converged, na.rm=TRUE)))
lg(sprintf("  전치 비율 범위      : %.1f%% ~ %.1f%%",
           100 * min(summary_df$prop_transposed, na.rm=TRUE),
           100 * max(summary_df$prop_transposed, na.rm=TRUE)))
flush_log()
cat(sprintf("\n로그 저장 완료: %s\n", LOG_PATH))
