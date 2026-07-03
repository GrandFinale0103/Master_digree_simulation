# =============================================================================
# 07_analysis_transposition.R
# 경계 모수 전치(boundary parameter transposition) 분석
#
# 분석 대상: 문항 1 (조작 문항)의 추정 경계 모수 b1 < b2 < b3 < b4 순서 유지 여부
#
# 출력 파일:
#   output/analysis/transposition_long.csv      — 반복별 long format 원자료
#   output/analysis/transposition_summary.csv   — 조건별 집계 요약
#   output/analysis/transposition_glm.txt       — GLM 분석 결과 (효과 검정 + 대비)
#   output/analysis/transposition_plot.png      — 조건별 전치 비율 시각화
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1단계: estimated_params CSV 파일 로딩 및 통합
# =============================================================================
est_files <- list.files(
  path       = "output/estimated_params",
  pattern    = "_est_params\\.csv$",
  full.names = TRUE
)

if (length(est_files) == 0) stop("추정 결과 파일이 없습니다. 먼저 시뮬레이션을 실행하세요.")

cat(sprintf("파일 %d개 로딩 중...\n", length(est_files)))

parse_fname <- function(path) {
  fname <- basename(path)
  cond  <- str_match(fname, "_cond([0-9]+)_")[, 2]
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
    cond_code = meta$cond_code,
    rep_id    = meta$rep_id,
    converged = isTRUE(item1$converged[1]),
    b1        = item1$b1[1],
    b2        = item1$b2[1],
    b3        = item1$b3[1],
    b4        = item1$b4[1],
    stringsAsFactors = FALSE
  )
}

raw_list <- lapply(est_files, read_est_file)
raw      <- do.call(rbind, raw_list[!sapply(raw_list, is.null)])
cat(sprintf("  총 %d행 로딩 완료\n", nrow(raw)))

# =============================================================================
# 2단계: 조건 코드 분해 (IV1–IV4)
# =============================================================================
# IV4 표기:
#   내부 코드(iv4_theta_dist): "pos_skew" / "normal" / "neg_skew" / "uniform"
#   표시 라벨(iv4_label)     : "정적편포(+0.8)" / "정규분포" / "부적편포(-0.8)" / "균등분포"

IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 2.0,  "2" = 1.5,  "3" = 1.0)
IV3_MAP <- c("1" = 0.0,  "2" = 1.5,  "3" = 3.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")
IV4_LABEL <- c(
  "pos_skew" = "정적편포(+0.8)",
  "normal"   = "정규분포",
  "neg_skew" = "부적편포(-0.8)",
  "uniform"  = "균등분포"
)

raw <- raw %>%
  mutate(
    d1 = substr(cond_code, 1, 1),
    d2 = substr(cond_code, 2, 2),
    d3 = substr(cond_code, 3, 3),
    d4 = substr(cond_code, 4, 4),
    iv1_sf4        = IV1_MAP[d1],
    iv2_b_interval = IV2_MAP[d2],
    iv3_b_mean     = IV3_MAP[d3],
    iv4_theta_dist = IV4_MAP[d4],
    iv4_label      = IV4_LABEL[IV4_MAP[d4]]
  ) %>%
  select(-d1, -d2, -d3, -d4)

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
         iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist, iv4_label,
         converged, b1, b2, b3, b4, transposed, which_transposed) %>%
  arrange(cond_code, rep_id)

write.csv(long_df, "output/analysis/transposition_long.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/transposition_long.csv (%d행)\n", nrow(long_df)))

# =============================================================================
# 5단계: 조건별 집계
# =============================================================================
summary_df <- long_df %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean,
           iv4_theta_dist, iv4_label) %>%
  summarise(
    n_total         = n(),
    n_converged     = sum(converged,   na.rm = TRUE),
    n_transposed    = sum(transposed,  na.rm = TRUE),
    prop_transposed = n_transposed / n_converged,
    n_12 = sum(grepl("12", which_transposed), na.rm = TRUE),
    n_23 = sum(grepl("23", which_transposed), na.rm = TRUE),
    n_34 = sum(grepl("34", which_transposed), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, "output/analysis/transposition_summary.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/transposition_summary.csv (%d행)\n",
            nrow(summary_df)))

# =============================================================================
# 6단계: GLM 분석
# =============================================================================
glm_data <- long_df %>%
  filter(converged == TRUE, !is.na(transposed)) %>%
  mutate(
    iv1 = factor(as.character(iv1_sf4),
                 levels = c("3", "3.33", "3.66"),
                 labels = c("sf4=3.00", "sf4=3.33", "sf4=3.66")),
    iv2 = factor(as.character(iv2_b_interval),
                 levels = c("2", "1.5", "1"),
                 labels = c("간격=2.0", "간격=1.5", "간격=1.0")),
    iv3 = factor(as.character(iv3_b_mean),
                 levels = c("0", "1.5", "3"),
                 labels = c("심각도=0", "심각도=1.5", "심각도=3.0")),
    iv4 = factor(iv4_theta_dist,
                 levels = c("normal", "pos_skew", "neg_skew", "uniform"),
                 labels = c("정규분포", "정적편포(+)", "부적편포(-)", "균등분포"))
  )

cat(sprintf("\nGLM 분석 대상: %d 반복 (수렴 성공 기준)\n", nrow(glm_data)))

# ── 완전 요인 모형: 주효과 + 2원(6) + 3원(4) + 4원(1) ────────────────────────
glm_full <- glm(
  transposed ~ iv1 * iv2 * iv3 * iv4,
  data   = glm_data,
  family = binomial(link = "logit")
)

# ── 편부분 pseudo-R² 계산 함수 ────────────────────────────────────────────────
# 각 항을 제거한 모형과 전체 모형의 이탈도 차이 / 절편만 모형의 이탈도
partial_pseudo_r2 <- function(full_model, term, data) {
  formula_reduced <- update(formula(full_model), paste(". ~ . -", term))
  fit_red <- tryCatch(
    glm(formula_reduced, data = data, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  if (is.null(fit_red)) return(NA_real_)
  dev_null <- full_model$null.deviance
  delta_dev <- deviance(fit_red) - deviance(full_model)
  delta_dev / dev_null
}

# ── Type III Wald 검정표 ───────────────────────────────────────────────────────
if (!requireNamespace("car", quietly = TRUE)) {
  install.packages("car", repos = "https://cran.rstudio.com/")
}
anova_tbl <- car::Anova(glm_full, type = 3)

# 효과크기 계산 (편부분 pseudo-R²)
terms_to_test <- rownames(anova_tbl)[rownames(anova_tbl) != "(Intercept)"]
cat("효과크기 계산 중 (항별 편부분 pseudo-R²)...\n")
partial_r2 <- sapply(terms_to_test, partial_pseudo_r2,
                     full_model = glm_full, data = glm_data)

# ── 대비 검정 (emmeans) ────────────────────────────────────────────────────────
if (!requireNamespace("emmeans", quietly = TRUE)) {
  install.packages("emmeans", repos = "https://cran.rstudio.com/")
}

emm_iv1 <- emmeans::emmeans(glm_full, ~ iv1, type = "response")
emm_iv2 <- emmeans::emmeans(glm_full, ~ iv2, type = "response")
emm_iv3 <- emmeans::emmeans(glm_full, ~ iv3, type = "response")
emm_iv4 <- emmeans::emmeans(glm_full, ~ iv4, type = "response")

contrast_iv1 <- emmeans::contrast(emm_iv1, method = "pairwise", adjust = "bonferroni")
contrast_iv2 <- emmeans::contrast(emm_iv2, method = "pairwise", adjust = "bonferroni")
contrast_iv3 <- emmeans::contrast(emm_iv3, method = "pairwise", adjust = "bonferroni")
contrast_iv4 <- emmeans::contrast(emm_iv4, method = "pairwise", adjust = "bonferroni")

# =============================================================================
# GLM 결과 저장
# =============================================================================
sink("output/analysis/transposition_glm.txt")

cat("================================================================\n")
cat("GLM 분석: 경계 모수 전치 여부 ~ 조건\n")
cat("family = binomial(logit) | 완전 요인 설계\n")
cat("================================================================\n\n")

cat("── 표본 크기 ──\n")
cat(sprintf("  전체: %d반복  수렴 성공: %d반복  전치 발생: %d반복 (%.1f%%)\n\n",
            nrow(long_df),
            sum(long_df$converged, na.rm=TRUE),
            sum(long_df$transposed, na.rm=TRUE),
            100 * mean(long_df$transposed, na.rm=TRUE)))

cat("── 모형 적합도 ──\n")
cat(sprintf("  Null deviance    : %.2f (df = %d)\n",
            glm_full$null.deviance, glm_full$df.null))
cat(sprintf("  Residual deviance: %.2f (df = %d)\n",
            deviance(glm_full), glm_full$df.residual))
cat(sprintf("  McFadden pseudo-R²: %.4f\n",
            1 - deviance(glm_full) / glm_full$null.deviance))
cat(sprintf("  AIC: %.2f\n\n", AIC(glm_full)))

cat("── Type III Wald 검정 + 편부분 pseudo-R² ──\n")
cat("  ※ 각 항을 제거했을 때 이탈도 증가분 / 절편 모형 이탈도\n\n")

# ANOVA 표에 pseudo-R² 열 추가
anova_out <- as.data.frame(anova_tbl)
anova_out$partial_pseudo_R2 <- NA_real_
anova_out[terms_to_test, "partial_pseudo_R2"] <- partial_r2

# 항 구분 추가
anova_out$effect_order <- NA_character_
for (nm in terms_to_test) {
  n_colon <- str_count(nm, ":")
  anova_out[nm, "effect_order"] <- switch(as.character(n_colon),
    "0" = "주효과",
    "1" = "2원 상호작용",
    "2" = "3원 상호작용",
    "3" = "4원 상호작용"
  )
}

print(anova_out[terms_to_test, ], digits = 4)
cat("\n")

cat("── 오즈비 (Odds Ratio) 및 95% CI — 계수 수준 ──\n")
or_ci <- tryCatch(
  exp(cbind(OR = coef(glm_full), confint(glm_full))),
  error = function(e) exp(cbind(OR = coef(glm_full),
                                 confint.default(glm_full)))
)
print(round(or_ci, 4))
cat("\n")

cat("────────────────────────────────────────────────────────────────\n")
cat("대비 검정 (Bonferroni 교정, 반응 확률 척도)\n")
cat("────────────────────────────────────────────────────────────────\n\n")

cat("── IV1 (채점함수 sf4) 주변 평균 ──\n")
print(emm_iv1)
cat("\n── IV1 쌍별 대비 ──\n")
print(contrast_iv1)

cat("\n── IV2 (경계모수 간격) 주변 평균 ──\n")
print(emm_iv2)
cat("\n── IV2 쌍별 대비 ──\n")
print(contrast_iv2)

cat("\n── IV3 (문항 심각도) 주변 평균 ──\n")
print(emm_iv3)
cat("\n── IV3 쌍별 대비 ──\n")
print(contrast_iv3)

cat("\n── IV4 (능력모수 분포) 주변 평균 ──\n")
print(emm_iv4)
cat("\n── IV4 쌍별 대비 ──\n")
print(contrast_iv4)

sink()
cat("저장 완료: output/analysis/transposition_glm.txt\n")

# =============================================================================
# 7단계: 시각화
# =============================================================================
plot_df <- summary_df %>%
  mutate(
    iv4_label = factor(iv4_label,
                       levels = c("정규분포", "정적편포(+0.8)",
                                  "부적편포(-0.8)", "균등분포"))
  )

p <- ggplot(plot_df,
            aes(x    = factor(iv3_b_mean),
                y    = prop_transposed,
                fill = iv4_label)) +
  geom_col(position = "dodge") +
  facet_grid(
    rows = vars(iv1_sf4),
    cols = vars(iv2_b_interval),
    labeller = labeller(
      iv1_sf4        = function(x) paste0("sf4 = ", x),
      iv2_b_interval = function(x) paste0("b 간격 = ", x)
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title = "조건별 경계 모수 전치 발생 비율",
    x     = "문항 심각도 — IV3: 경계모수 평균 (0 / 1.5 / 3.0)",
    y     = "전치 발생 비율",
    fill  = "능력모수 분포 (IV4)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/analysis/transposition_plot.png", p,
       width = 10, height = 8, dpi = 150)
cat("저장 완료: output/analysis/transposition_plot.png\n")

# =============================================================================
# 최종 요약
# =============================================================================
cat("\n================================================================\n")
cat("전체 분석 완료\n")
cat("================================================================\n")
cat(sprintf("  전체 반복 수        : %d\n", nrow(long_df)))
cat(sprintf("  수렴 성공           : %d (%.1f%%)\n",
            sum(long_df$converged, na.rm=TRUE),
            100 * mean(long_df$converged, na.rm=TRUE)))
cat(sprintf("  전치 발생 (전체)    : %d (%.1f%%)\n",
            sum(long_df$transposed, na.rm=TRUE),
            100 * mean(long_df$transposed, na.rm=TRUE)))
cat(sprintf("  전치 비율 범위      : %.1f%% ~ %.1f%%\n",
            100 * min(summary_df$prop_transposed, na.rm=TRUE),
            100 * max(summary_df$prop_transposed, na.rm=TRUE)))
cat("\n출력 파일:\n")
cat("  output/analysis/transposition_long.csv\n")
cat("  output/analysis/transposition_summary.csv\n")
cat("  output/analysis/transposition_glm.txt\n")
cat("  output/analysis/transposition_plot.png\n")
