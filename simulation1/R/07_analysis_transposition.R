# =============================================================================
# 07_analysis_transposition.R
# 경계 모수 전치(boundary parameter transposition) 분석
#
# 분석 대상: 문항 1 (조작 문항)의 추정 경계 모수 b1 < b2 < b3 < b4 순서 유지 여부
#
# 출력 파일:
#   output/analysis/transposition_long.csv      — 반복별 long format 원자료
#   output/analysis/transposition_summary.csv   — 조건별 집계 요약
#   output/analysis/transposition_glm.txt       — GLM 분석 결과
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
# 파일명 형식: YYYYMMDD_cond[C]_rep[R]_seed[S]_est_params.csv

est_files <- list.files(
  path       = "output/estimated_params",
  pattern    = "_est_params\\.csv$",
  full.names = TRUE
)

if (length(est_files) == 0) stop("추정 결과 파일이 없습니다. 먼저 시뮬레이션을 실행하세요.")

cat(sprintf("파일 %d개 로딩 중...\n", length(est_files)))

# 파일명에서 조건 코드와 반복 번호 파싱
parse_fname <- function(path) {
  fname <- basename(path)
  cond  <- str_match(fname, "_cond([0-9]+)_")[, 2]
  rep   <- as.integer(str_match(fname, "_rep([0-9]+)_")[, 2])
  list(cond_code = cond, rep_id = rep)
}

# 각 파일에서 문항 1 행만 추출
read_est_file <- function(path) {
  meta <- parse_fname(path)
  df   <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (is.null(df)) return(NULL)

  # 추정 실패 파일 처리 (success=FALSE 열이 있는 경우)
  if ("success" %in% names(df) && isFALSE(df$success[1])) {
    return(data.frame(
      cond_code       = meta$cond_code,
      rep_id          = meta$rep_id,
      converged       = FALSE,
      b1 = NA_real_, b2 = NA_real_, b3 = NA_real_, b4 = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # 문항 1 행 추출 (item 열 기준)
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
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 2.0,  "2" = 1.5,  "3" = 1.0)
IV3_MAP <- c("1" = 0.0,  "2" = 1.5,  "3" = 3.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")

raw <- raw %>%
  mutate(
    d1 = substr(cond_code, 1, 1),
    d2 = substr(cond_code, 2, 2),
    d3 = substr(cond_code, 3, 3),
    d4 = substr(cond_code, 4, 4),
    iv1_sf4       = IV1_MAP[d1],
    iv2_b_interval= IV2_MAP[d2],
    iv3_b_mean    = IV3_MAP[d3],
    iv4_theta_dist= IV4_MAP[d4]
  ) %>%
  select(-d1, -d2, -d3, -d4)

# =============================================================================
# 3단계: 경계 모수 전치 탐지
# =============================================================================
# 전치 기준: b_k >= b_{k+1} (같거나 역전된 경우 모두 전치로 간주)
# which_transposed: 위반된 인접 쌍의 인덱스를 연결한 문자열
#   예) b2 >= b3 → "23", b1 >= b2 & b3 >= b4 → "1234"

detect_transposition <- function(b1, b2, b3, b4) {
  if (any(is.na(c(b1, b2, b3, b4)))) {
    return(list(transposed = NA_integer_, which_transposed = NA_character_))
  }
  pairs     <- list(c(1,2), c(2,3), c(3,4))
  bvals     <- c(b1, b2, b3, b4)
  violated  <- sapply(pairs, function(p) bvals[p[1]] >= bvals[p[2]])
  which_str <- paste0(
    sapply(which(violated), function(i) paste0(pairs[[i]], collapse = "")),
    collapse = ""
  )
  list(
    transposed      = as.integer(any(violated)),
    which_transposed = if (any(violated)) which_str else ""
  )
}

trans_results <- mapply(
  detect_transposition,
  raw$b1, raw$b2, raw$b3, raw$b4,
  SIMPLIFY = FALSE
)

raw$transposed       <- sapply(trans_results, `[[`, "transposed")
raw$which_transposed <- sapply(trans_results, `[[`, "which_transposed")

# =============================================================================
# 4단계: Long format CSV 저장
# =============================================================================
long_df <- raw %>%
  select(
    cond_code, rep_id,
    iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
    converged,
    b1, b2, b3, b4,
    transposed, which_transposed
  ) %>%
  arrange(cond_code, rep_id)

write.csv(long_df, "output/analysis/transposition_long.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/transposition_long.csv (%d행)\n", nrow(long_df)))

# =============================================================================
# 5단계: 조건별 집계
# =============================================================================
summary_df <- long_df %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist) %>%
  summarise(
    n_total       = n(),
    n_converged   = sum(converged, na.rm = TRUE),
    n_transposed  = sum(transposed, na.rm = TRUE),
    prop_transposed = n_transposed / n_converged,
    n_12 = sum(grepl("12", which_transposed), na.rm = TRUE),
    n_23 = sum(grepl("23", which_transposed), na.rm = TRUE),
    n_34 = sum(grepl("34", which_transposed), na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, "output/analysis/transposition_summary.csv", row.names = FALSE)
cat(sprintf("저장 완료: output/analysis/transposition_summary.csv (%d행)\n", nrow(summary_df)))

# =============================================================================
# 6단계: GLM 분석 (binomial logistic regression)
# =============================================================================
# 수렴 성공 + 전치 여부가 NA가 아닌 반복만 사용
glm_data <- long_df %>%
  filter(converged == TRUE, !is.na(transposed)) %>%
  mutate(
    iv1 = factor(iv1_sf4,        levels = c(3.00, 3.33, 3.66)),
    iv2 = factor(iv2_b_interval, levels = c(2.0, 1.5, 1.0)),
    iv3 = factor(iv3_b_mean,     levels = c(0.0, 1.5, 3.0)),
    iv4 = factor(iv4_theta_dist, levels = c("normal", "pos_skew", "neg_skew", "uniform"))
  )

cat(sprintf("\nGLM 분석 대상: %d 반복 (수렴 성공 기준)\n", nrow(glm_data)))

# 주효과 + 2차 상호작용 모형
glm_fit <- glm(
  transposed ~ iv1 + iv2 + iv3 + iv4 +
    iv1:iv2 + iv1:iv3 + iv1:iv4 +
    iv2:iv3 + iv2:iv4 +
    iv3:iv4,
  data   = glm_data,
  family = binomial(link = "logit")
)

# 결과 저장
sink("output/analysis/transposition_glm.txt")
cat("========================================\n")
cat("GLM 분석: 경계 모수 전치 여부 ~ 조건\n")
cat("family = binomial(logit)\n")
cat("========================================\n\n")
cat("── 모형 요약 ──\n")
print(summary(glm_fit))
cat("\n── 분산 분석 (Wald chi-square, type III) ──\n")
if (requireNamespace("car", quietly = TRUE)) {
  print(car::Anova(glm_fit, type = 3))
} else {
  print(anova(glm_fit, test = "Chisq"))
}
cat("\n── 오즈비 (Odds Ratio) 및 95% CI ──\n")
or_ci <- exp(cbind(OR = coef(glm_fit), confint(glm_fit)))
print(round(or_ci, 4))
sink()

cat("저장 완료: output/analysis/transposition_glm.txt\n")

# =============================================================================
# 7단계: 시각화
# =============================================================================
# 조건별 전치 비율 barplot (IV1 × IV2 × IV3 패널, IV4 색상)
p <- ggplot(summary_df,
            aes(x    = factor(iv3_b_mean),
                y    = prop_transposed,
                fill = iv4_theta_dist)) +
  geom_col(position = "dodge") +
  facet_grid(iv1_sf4 ~ iv2_b_interval,
             labeller = labeller(
               iv1_sf4        = function(x) paste0("sf4=", x),
               iv2_b_interval = function(x) paste0("간격=", x)
             )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(
    title = "조건별 경계 모수 전치 발생 비율",
    x     = "문항 심각도 (IV3: 경계모수 평균)",
    y     = "전치 발생 비율",
    fill  = "능력모수 분포 (IV4)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave("output/analysis/transposition_plot.png", p,
       width = 10, height = 8, dpi = 150)
cat("저장 완료: output/analysis/transposition_plot.png\n")

# =============================================================================
# 요약 출력
# =============================================================================
cat("\n========================================\n")
cat("전체 분석 완료\n")
cat("========================================\n")
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
