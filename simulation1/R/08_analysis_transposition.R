# =============================================================================
# 08_analysis_transposition.R
# 경계 모수 전치(boundary parameter transposition) 분석
#
# [파일 검증]
#   - 기대 조건 코드(IV 수준에서 자동 산출) 대비 실제 파일 존재 여부를 검사
#   - 누락 조건, 반복 횟수 미달 조건을 로그에 기록
#   - 존재하는 파일만으로 유연하게 집계
#
# 분석 대상: 문항 1 (조작 문항)의 추정 경계 모수 b1 < b2 < b3 < b4 순서 유지 여부
#
# 출력 파일:
#   output/analysis/08_log.txt                  — 실행 로그
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

# ── 로그 시스템 ───────────────────────────────────────────────────────────────
LOG_PATH  <- "output/analysis/08_log.txt"
log_lines <- character(0)
lg <- function(...) {
  msg <- paste0(...); cat(msg, "\n"); log_lines <<- c(log_lines, msg)
}
lg_section <- function(title) lg(sprintf("\n[%s]  %s", Sys.time(), title))
flush_log  <- function() writeLines(log_lines, LOG_PATH)

# ── 조건 코드 → 실제 값 매핑 (수준값·수준 수 변경 시 이 블록만 수정) ─────────
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)   # 경계모수 간격
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)   # 문항 심각도
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")
IV4_LABEL <- c(
  "pos_skew" = "정적편포(+0.8)",
  "normal"   = "정규분포",
  "neg_skew" = "부적편포(-0.8)",
  "uniform"  = "균등분포"
)

# 기대 조건 코드 전체 목록 (IV 수준 수에서 자동 생성)
EXPECTED_CONDS <- character(0)
for (d1 in seq_along(IV1_MAP)) for (d2 in seq_along(IV2_MAP))
  for (d3 in seq_along(IV3_MAP)) for (d4 in seq_along(IV4_MAP))
    EXPECTED_CONDS <- c(EXPECTED_CONDS, paste0(d1, d2, d3, d4))
N_CONDS <- length(EXPECTED_CONDS)  # 현재 설계: 108

# ── 파일 검증 함수 ────────────────────────────────────────────────────────────
validate_est_files <- function(files) {
  lg_section("파일 검증 (est_params)")
  lg(sprintf("  발견 파일 수: %d", length(files)))

  parse_cond <- function(p) {
    m <- regmatches(p, regexpr("_cond([0-9]+)_", p))
    if (length(m) == 0) return(NA_character_)
    gsub("_cond|_", "", m)
  }

  found_conds <- sapply(files, parse_cond)
  reps_tbl    <- sort(table(found_conds[!is.na(found_conds)]))

  missing_conds <- setdiff(EXPECTED_CONDS, names(reps_tbl))
  extra_conds   <- setdiff(names(reps_tbl), EXPECTED_CONDS)

  if (length(missing_conds) > 0) {
    lg(sprintf("  [경고] 누락 조건 %d개: %s",
               length(missing_conds), paste(missing_conds, collapse = ", ")))
  } else {
    lg(sprintf("  누락 조건 없음 (%d개 전체 존재)", N_CONDS))
  }

  if (length(extra_conds) > 0) {
    lg(sprintf("  [경고] 예상 외 조건 %d개: %s",
               length(extra_conds), paste(extra_conds, collapse = ", ")))
  }

  if (length(reps_tbl) > 0) {
    max_reps    <- max(reps_tbl)
    under_conds <- reps_tbl[reps_tbl < max_reps]
    lg(sprintf("  최대 반복 횟수: %d", max_reps))
    if (length(under_conds) > 0) {
      lg(sprintf("  [경고] 반복 미달 조건 %d개 (기준: %d회):",
                 length(under_conds), max_reps))
      for (nm in names(under_conds)) {
        lg(sprintf("    cond=%s: %d회", nm, under_conds[[nm]]))
      }
    } else {
      lg(sprintf("  모든 조건이 동일 반복 횟수(%d회) 충족", max_reps))
    }
  }
  invisible(NULL)
}

# =============================================================================
# 1단계: estimated_params CSV 파일 로딩 및 통합
# =============================================================================
lg_section("1단계: est_params 파일 로딩")

est_files <- list.files(
  path       = "output/estimated_params",
  pattern    = "_est_params\\.csv$",
  full.names = TRUE
)

if (length(est_files) == 0) {
  flush_log()
  stop("추정 결과 파일이 없습니다. 먼저 시뮬레이션을 실행하세요.")
}

validate_est_files(est_files)

lg(sprintf("파일 %d개 로딩 중...", length(est_files)))

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
    converged = as.logical(item1$converged[1]),   # isTRUE("TRUE")=FALSE 방지
    b1        = item1$b1[1],
    b2        = item1$b2[1],
    b3        = item1$b3[1],
    b4        = item1$b4[1],
    stringsAsFactors = FALSE
  )
}

raw_list <- lapply(est_files, read_est_file)
raw      <- do.call(rbind, raw_list[!sapply(raw_list, is.null)])
lg(sprintf("  총 %d행 로딩 완료", nrow(raw)))

# =============================================================================
# 2단계: 조건 코드 분해 (IV1–IV4)
# =============================================================================
lg_section("2단계: 조건 코드 분해")
# IV4 표기:
#   내부 코드(iv4_theta_dist): "pos_skew" / "normal" / "neg_skew" / "uniform"
#   표시 라벨(iv4_label)     : "정적편포(+0.8)" / "정규분포" / "부적편포(-0.8)" / "균등분포"

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
lg(sprintf("저장 완료: output/analysis/transposition_long.csv (%d행)", nrow(long_df)))

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
lg(sprintf("저장 완료: output/analysis/transposition_summary.csv (%d행)", nrow(summary_df)))

# =============================================================================
# 6단계: GLM 분석
# =============================================================================

# ── helper 1: 수치형 IV → factor (데이터 관찰값 기반, 설계 변경에 자동 적응) ─
make_iv_factor <- function(x, fmt = "%.2f", decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}

# ── helper 2: 문자형 IV → factor (정의된 수준 중 데이터에 존재하는 것만 사용) ─
make_char_factor <- function(x, ordered_levels, labels = NULL) {
  obs_lv <- ordered_levels[ordered_levels %in% unique(na.omit(x))]
  obs_lb <- if (!is.null(labels)) labels[match(obs_lv, ordered_levels)] else obs_lv
  factor(x, levels = obs_lv, labels = obs_lb)
}

# ── IV 설정 테이블: 여기만 수정하면 factor 변환·검정이 자동 적용됨 ─────────
IV_CONFIG <- list(
  iv1 = list(col = "iv1_sf4",        type = "numeric", fmt = "%.2f", decreasing = FALSE),
  iv2 = list(col = "iv2_b_interval", type = "numeric", fmt = "%.1f", decreasing = TRUE),
  iv3 = list(col = "iv3_b_mean",     type = "numeric", fmt = "%.1f", decreasing = FALSE),
  iv4 = list(col = "iv4_theta_dist", type = "char",
             levels = c("normal", "pos_skew", "neg_skew", "uniform"),
             labels = c("정규분포", "정적편포(+)", "부적편포(-)", "균등분포"))
)

# ── 데이터 준비 ───────────────────────────────────────────────────────────────
glm_data <- long_df %>%
  filter(converged == TRUE, !is.na(transposed)) %>%
  mutate(transposed = as.numeric(transposed))   # binomial GLM: 반드시 numeric 0/1

# ── 자동 factor 변환 (IV_CONFIG 루프) ────────────────────────────────────────
for (iv_name in names(IV_CONFIG)) {
  cfg <- IV_CONFIG[[iv_name]]
  glm_data[[iv_name]] <-
    if (cfg$type == "numeric") {
      make_iv_factor(glm_data[[cfg$col]], fmt = cfg$fmt, decreasing = cfg$decreasing)
    } else {
      make_char_factor(glm_data[[cfg$col]], cfg$levels, cfg$labels)
    }
}

iv_vars <- names(IV_CONFIG)

# ── 진단 1: 기본 데이터 확인 ─────────────────────────────────────────────────
lg_section("6단계: GLM 분석")
lg(sprintf("GLM 분석 대상: %d 반복 (수렴 성공 기준)", nrow(glm_data)))
lg(sprintf("  [진단] converged 분포: TRUE=%d, FALSE=%d, NA=%d",
            sum(long_df$converged == TRUE,  na.rm = TRUE),
            sum(long_df$converged == FALSE, na.rm = TRUE),
            sum(is.na(long_df$converged))))

if (nrow(glm_data) == 0) {
  flush_log()
  stop(paste(
    "수렴 성공 데이터가 없습니다.",
    "→ 추정 결과 파일의 converged 컬럼 또는 시뮬레이션 실행 여부를 확인하세요."
  ))
}
lg(sprintf("  [진단] transposed 분포: 0=%d건, 1=%d건, NA=%d건",
            sum(glm_data$transposed == 0, na.rm = TRUE),
            sum(glm_data$transposed == 1, na.rm = TRUE),
            sum(is.na(glm_data$transposed))))

if (isTRUE(var(glm_data$transposed, na.rm = TRUE) == 0)) {
  flush_log()
  stop("transposed 분산이 0입니다 (모두 동일) — GLM 적합 불가.")
}

lg("  IV별 실제 관측 수준:")
for (v in iv_vars) {
  lg(sprintf("    %s: %d개 수준 (%s)",
              v, nlevels(glm_data[[v]]),
              paste(levels(glm_data[[v]]), collapse = ", ")))
}

# ── 진단 2: 완전 분리(complete separation) 탐지 ──────────────────────────────
#
# [핵심 진단] Wald 검정(car::Anova 기본값)이 χ²=0, p=1을 반환하는 주요 원인:
#   완전 분리 — 일부 셀의 전치 발생률이 정확히 0% 또는 100%이면
#   logistic regression 계수가 ±∞로 발산하고 표준오차 → ∞가 됩니다.
#   Wald χ² = (추정치/SE)² → 0 / ∞ = 0, p → 1.
#   → 해결: 우도비 검정(LRT)은 계수 추정값이 아닌 이탈도 차이를 사용하므로 안정적.
#
cell_check <- glm_data %>%
  group_by(across(all_of(iv_vars))) %>%
  summarise(n = n(), prop_trans = mean(transposed, na.rm = TRUE), .groups = "drop")

n_sep_zero <- sum(cell_check$prop_trans == 0,   na.rm = TRUE)
n_sep_one  <- sum(cell_check$prop_trans == 1,   na.rm = TRUE)
lg(sprintf("  [완전 분리 진단] 전치 0%% 셀: %d개, 100%% 셀: %d개 (총 %d셀)",
            n_sep_zero, n_sep_one, nrow(cell_check)))
if (n_sep_zero + n_sep_one > 0) {
  lg("  ※ 완전 분리 탐지 — Wald χ²=0/p=1 오류 원인 확인됨")
  lg("     우도비 검정(LRT)으로 대체합니다.")
} else {
  lg("  완전 분리 없음 — LRT로 진행합니다 (Wald 대비 안정성 우선).")
}

# ── GLM 모형 적합 ─────────────────────────────────────────────────────────────
multi_ivs  <- iv_vars[sapply(iv_vars, function(v) nlevels(glm_data[[v]]) >= 2)]
single_ivs <- setdiff(iv_vars, multi_ivs)

if (length(single_ivs) > 0) {
  lg(sprintf("  ※ 단일 수준 요인 — 모형 제외: %s", paste(single_ivs, collapse = ", ")))
}
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

# ── 우도비 검정 (LRT) — Wald 대체 ────────────────────────────────────────────
# drop1()은 각 항을 제거한 모형과 전체 모형의 이탈도 차이(LRT 통계량)를 계산함.
# Type III 원리(다른 항은 유지, 해당 항만 제거)와 동일하며, 완전 분리에 강건함.
anova_tbl <- drop1(glm_full, scope = ~., test = "Chisq")

# LRT 통계량 컬럼명 탐지 (R 버전별 차이 대응)
lrt_col  <- intersect(c("LRT", "Chisq"), colnames(anova_tbl))[1]
pval_col <- grep("^Pr", colnames(anova_tbl), value = TRUE)[1]

# 편부분 pseudo-R²: LRT 통계량 / 절편 모형 이탈도
# (= 항 제거 시 이탈도 증가분 / 기저 이탈도)
lrt_terms  <- rownames(anova_tbl)[rownames(anova_tbl) != "<none>"]
partial_r2 <- setNames(
  anova_tbl[lrt_terms, lrt_col] / glm_full$null.deviance,
  lrt_terms
)

# ── 대비 검정 (emmeans) ────────────────────────────────────────────────────────
if (!requireNamespace("emmeans", quietly = TRUE)) {
  install.packages("emmeans", repos = "https://cran.rstudio.com/")
}

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

cat("── Type III 우도비 검정(LRT) + 편부분 pseudo-R² ──\n")
cat("  ※ Wald 검정은 완전 분리 시 χ²=0/p=1 오류 → drop1() LRT로 대체\n")
cat("  ※ 편부분 pseudo-R² = LRT 통계량 / 절편 모형 이탈도\n\n")

# drop1 출력에 효과 구분·pseudo-R² 열 추가
anova_out <- as.data.frame(anova_tbl[lrt_terms, , drop = FALSE])
anova_out$partial_pseudo_R2 <- partial_r2[lrt_terms]
anova_out$effect_order <- sapply(lrt_terms, function(nm) {
  nc <- str_count(nm, ":")
  switch(as.character(nc),
    "0" = "주효과", "1" = "2원 상호작용",
    "2" = "3원 상호작용", "3" = "4원 상호작용", nm)
})

print(anova_out, digits = 4)
cat("\n")

cat("── 완전 분리 진단 ──\n")
cat(sprintf("  전치 0%% 셀: %d개,  100%% 셀: %d개  (전체 %d셀)\n",
            n_sep_zero, n_sep_one, nrow(cell_check)))
if (n_sep_zero + n_sep_one > 0) {
  cat("  → OR/SE는 해당 셀에서 무한대로 발산할 수 있으므로 해석에 주의\n")
  cat("     (LRT 통계량은 이 문제에 영향받지 않음)\n")
}
cat("\n")

cat("── 오즈비 (Odds Ratio) 및 95% CI — 계수 수준 (Wald) ──\n")
cat("  ※ confint.default() = Wald (표준오차 기반, 즉시 계산)\n")
cat("  ※ 완전 분리 셀의 OR/CI는 극단값(Inf)이 될 수 있음\n\n")
or_ci <- exp(cbind(OR = coef(glm_full), confint.default(glm_full)))
print(round(or_ci, 4))
cat("\n")

cat("────────────────────────────────────────────────────────────────\n")
cat("대비 검정 (Bonferroni 교정, 반응 확률 척도)\n")
cat("────────────────────────────────────────────────────────────────\n\n")

iv_labels <- c(iv1 = "IV1 (채점함수 sf4)",
               iv2 = "IV2 (경계모수 간격)",
               iv3 = "IV3 (문항 심각도)",
               iv4 = "IV4 (능력모수 분포)")
for (v in multi_ivs) {
  cat(sprintf("\n── %s 주변 평균 ──\n", iv_labels[v]))
  print(emm_list[[v]])
  cat(sprintf("\n── %s 쌍별 대비 ──\n", iv_labels[v]))
  print(contrast_list[[v]])
}
for (v in single_ivs) {
  cat(sprintf("\n── %s: 수준이 1개뿐이어서 대비 생략 ──\n", iv_labels[v]))
}

sink()
lg("저장 완료: output/analysis/transposition_glm.txt")

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
lg("저장 완료: output/analysis/transposition_plot.png")

# =============================================================================
# 최종 요약
# =============================================================================
lg_section("최종 요약")
lg("================================================================")
lg("전체 분석 완료")
lg("================================================================")
lg(sprintf("  전체 반복 수        : %d", nrow(long_df)))
lg(sprintf("  수렴 성공           : %d (%.1f%%)",
            sum(long_df$converged, na.rm=TRUE),
            100 * mean(long_df$converged, na.rm=TRUE)))
lg(sprintf("  전치 발생 (전체)    : %d (%.1f%%)",
            sum(long_df$transposed, na.rm=TRUE),
            100 * mean(long_df$transposed, na.rm=TRUE)))
lg(sprintf("  전치 비율 범위      : %.1f%% ~ %.1f%%",
            100 * min(summary_df$prop_transposed, na.rm=TRUE),
            100 * max(summary_df$prop_transposed, na.rm=TRUE)))
lg("\n출력 파일:")
lg("  output/analysis/08_log.txt")
lg("  output/analysis/transposition_long.csv")
lg("  output/analysis/transposition_summary.csv")
lg("  output/analysis/transposition_glm.txt")
lg("  output/analysis/transposition_plot.png")

flush_log()
