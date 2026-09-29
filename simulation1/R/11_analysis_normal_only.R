# =============================================================================
# 11_analysis_normal_only.R
# 능력모수 분포(IV4)를 "정규분포"로 고정한 조건만 따로 분석
#
# 대상: 조건 코드 4번째 자리 = 2 (normal) → IV1(3) × IV2(3) × IV3(3) = 27 조건
# 독립변수: IV1(sf4) * IV2(경계모수 간격) * IV3(경계모수 평균)   ※ IV4 제외
#
# 기존 07~10 결과와 섞이지 않도록 모든 결과를 별도 폴더에 저장:
#   output/analysis_normal/
#     11_log.txt                    — 실행 로그
#     transposition_long.csv        — 반복별 전치 여부
#     transposition_summary.csv     — 조건별 전치 비율
#     transposition_glm.txt         — 전치 여부 로지스틱 GLM (LRT)
#     transposition_plot.png        — 조건별 전치 비율 그래프
#     rmse_bias_long.csv            — 반복별 오차 원자료
#     rmse_bias_summary.csv         — 조건 × 모수별 Bias / RMSE
#     rmse_bias_plot_b_bias.png     — b1~b4 Bias 그래프
#     rmse_bias_plot_b_rmse.png     — b1~b4 RMSE 그래프
#     bias_rmse_anova.txt           — b3/b4 Bias·RMSE Type III ANOVA (4개 표)
#     bias_ttest_by_cond.csv        — 조건별 Bias t검정 (FDR 보정)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})
if (!requireNamespace("car", quietly = TRUE))
  install.packages("car", repos = "https://cran.rstudio.com/")

OUT_DIR <- "output/analysis_normal"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

# ── 로그 ─────────────────────────────────────────────────────────────────────
LOG_PATH  <- out("11_log.txt")
log_lines <- character(0)
lg <- function(...) {
  msg <- paste0(...); cat(msg, "\n"); log_lines <<- c(log_lines, msg)
}
lg_section <- function(title) lg(sprintf("\n[%s]  %s", Sys.time(), title))
flush_log  <- function() writeLines(log_lines, LOG_PATH)

# ── 설계 상수 ────────────────────────────────────────────────────────────────
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
NORMAL_DIGIT <- "2"   # IV4: 1=pos_skew, 2=normal, 3=neg_skew, 4=uniform
DISCRIM <- 1

EXPECTED_CONDS <- as.vector(outer(
  outer(names(IV1_MAP), names(IV2_MAP), paste0), names(IV3_MAP), paste0))
EXPECTED_CONDS <- sort(paste0(EXPECTED_CONDS, NORMAL_DIGIT))   # 27개

lg("================================================================")
lg(sprintf("11_analysis_normal_only.R 실행 시작: %s", Sys.time()))
lg("분석 대상: IV4 = 정규분포 조건만 (IV4 고정)")
lg("================================================================")

# =============================================================================
# 1단계: 정규분포 조건 파일만 로딩 (문항 1)
# =============================================================================
lg_section("1단계: est_params 로딩 (정규분포 조건)")

all_files <- list.files("output/estimated_params",
                        pattern = "_est_params\\.csv$", full.names = TRUE)
file_cond <- str_match(basename(all_files), "_cond([0-9]+)_")[, 2]
est_files <- all_files[!is.na(file_cond) & substr(file_cond, 4, 4) == NORMAL_DIGIT]

if (length(est_files) == 0) {
  flush_log()
  stop("정규분포 조건의 추정 결과 파일이 없습니다.")
}

found <- table(str_match(basename(est_files), "_cond([0-9]+)_")[, 2])
missing <- setdiff(EXPECTED_CONDS, names(found))
lg(sprintf("  파일 수: %d  |  발견 조건: %d / %d",
           length(est_files), length(found), length(EXPECTED_CONDS)))
if (length(missing) > 0)
  lg(sprintf("  [경고] 누락 조건: %s", paste(missing, collapse = ", ")))
if (length(unique(found)) > 1)
  lg(sprintf("  [경고] 조건별 반복 수 불일치: %d ~ %d", min(found), max(found)))

read_item1_est <- function(path) {
  cond <- str_match(basename(path), "_cond([0-9]+)_")[, 2]
  rep  <- as.integer(str_match(basename(path), "_rep([0-9]+)_")[, 2])
  df   <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  if ("success" %in% names(df) && isFALSE(df$success[1])) {
    return(data.frame(cond_code = cond, rep_id = rep, converged = FALSE,
                      a_est = NA_real_, b1_est = NA_real_, b2_est = NA_real_,
                      b3_est = NA_real_, b4_est = NA_real_))
  }
  item1 <- df[df$item %in% c("Item1", "item1", 1), ]
  if (nrow(item1) == 0) item1 <- df[1, ]
  data.frame(cond_code = cond, rep_id = rep,
             converged = as.logical(item1$converged[1]),
             a_est  = item1$a[1],
             b1_est = item1$b1[1], b2_est = item1$b2[1],
             b3_est = item1$b3[1], b4_est = item1$b4[1])
}

est_list <- lapply(est_files, read_item1_est)
raw <- do.call(rbind, est_list[!sapply(est_list, is.null)])
lg(sprintf("  총 %d행 로딩 완료", nrow(raw)))

# =============================================================================
# 2단계: 조건 분해, 참값·오차·전치 계산
# =============================================================================
raw <- raw %>%
  mutate(
    iv1_sf4        = unname(IV1_MAP[substr(cond_code, 1, 1)]),
    iv2_b_interval = unname(IV2_MAP[substr(cond_code, 2, 2)]),
    iv3_b_mean     = unname(IV3_MAP[substr(cond_code, 3, 3)]),
    a_true  = DISCRIM,
    b1_true = iv3_b_mean + iv2_b_interval * (-1.5),
    b2_true = iv3_b_mean + iv2_b_interval * (-0.5),
    b3_true = iv3_b_mean + iv2_b_interval *   0.5,
    b4_true = iv3_b_mean + iv2_b_interval *   1.5,
    err_a  = a_est  - a_true,
    err_b1 = b1_est - b1_true,
    err_b2 = b2_est - b2_true,
    err_b3 = b3_est - b3_true,
    err_b4 = b4_est - b4_true,
    t12 = as.integer(b1_est >= b2_est),
    t23 = as.integer(b2_est >= b3_est),
    t34 = as.integer(b3_est >= b4_est),
    transposed = as.integer(t12 | t23 | t34)
  ) %>%
  arrange(cond_code, rep_id)

conv_df <- raw %>% filter(converged == TRUE)
lg(sprintf("  수렴 성공: %d / %d", nrow(conv_df), nrow(raw)))

make_iv_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(
    f_iv1 = make_iv_factor(iv1_sf4,        "%.2f"),
    f_iv2 = make_iv_factor(iv2_b_interval, "%.1f"),
    f_iv3 = make_iv_factor(iv3_b_mean,     "%.1f")
  )
}

# =============================================================================
# 3단계: 경계모수 전치 분석
# =============================================================================
lg_section("3단계: 전치 분석")

trans_long <- raw %>%
  select(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean, converged,
         b1 = b1_est, b2 = b2_est, b3 = b3_est, b4 = b4_est,
         transposed, t12, t23, t34)
write.csv(trans_long, out("transposition_long.csv"), row.names = FALSE)

trans_summary <- trans_long %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
  summarise(
    n_total         = n(),
    n_converged     = sum(converged, na.rm = TRUE),
    n_transposed    = sum(transposed[converged %in% TRUE], na.rm = TRUE),
    prop_transposed = n_transposed / n_converged,
    n_12 = sum(t12[converged %in% TRUE], na.rm = TRUE),
    n_23 = sum(t23[converged %in% TRUE], na.rm = TRUE),
    n_34 = sum(t34[converged %in% TRUE], na.rm = TRUE),
    .groups = "drop"
  )
write.csv(trans_summary, out("transposition_summary.csv"), row.names = FALSE)
lg(sprintf("  저장 완료: transposition_long.csv / transposition_summary.csv (%d 조건)",
           nrow(trans_summary)))

# ── 로지스틱 GLM + LRT (완전 분리에 강건) ────────────────────────────────────
glm_data <- conv_df %>% filter(!is.na(transposed)) %>% add_factors()

sink(out("transposition_glm.txt"))
tryCatch({
  cat("================================================================\n")
  cat("GLM: 전치 여부 ~ IV1 * IV2 * IV3  (IV4 = 정규분포 고정)\n")
  cat("family = binomial(logit)  |  검정: drop1() 우도비 검정(LRT)\n")
  cat("================================================================\n\n")
  cat(sprintf("수렴 성공 반복: %d  |  전치 발생: %d (%.1f%%)\n\n",
              nrow(glm_data), sum(glm_data$transposed),
              100 * mean(glm_data$transposed)))
  if (var(glm_data$transposed) == 0) {
    cat("전치 여부 분산이 0 — GLM 적합 불가\n")
  } else {
    glm_full <- suppressWarnings(
      glm(transposed ~ f_iv1 * f_iv2 * f_iv3, data = glm_data,
          family = binomial(link = "logit")))
    lrt <- drop1(glm_full, scope = ~., test = "Chisq")
    lrt_col <- intersect(c("LRT", "Chisq"), colnames(lrt))[1]
    terms   <- rownames(lrt)[rownames(lrt) != "<none>"]
    lrt_out <- as.data.frame(lrt[terms, , drop = FALSE])
    lrt_out$partial_pseudo_R2 <- lrt_out[[lrt_col]] / glm_full$null.deviance
    cat(sprintf("McFadden pseudo-R²: %.4f\n\n",
                1 - deviance(glm_full) / glm_full$null.deviance))
    print(lrt_out, digits = 4)
  }
}, error = function(e) {
  cat("\n[오류] GLM 중 에러:", conditionMessage(e), "\n")
}, finally = sink())
lg("  저장 완료: transposition_glm.txt")

p_trans <- ggplot(add_factors(trans_summary),
                  aes(x = f_iv3, y = prop_transposed)) +
  geom_col(fill = "steelblue") +
  facet_grid(rows = vars(f_iv1), cols = vars(f_iv2),
             labeller = labeller(f_iv1 = function(x) paste0("s4 = ", x),
                                 f_iv2 = function(x) paste0("Interval = ", x))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "b-parameter Mean (IV3)", y = "Reversal Rate") +
  theme_bw(base_size = 14)
ggsave(out("transposition_plot.png"), p_trans, width = 10, height = 8, dpi = 150)
lg("  저장 완료: transposition_plot.png")

# =============================================================================
# 4단계: Bias / RMSE
# =============================================================================
lg_section("4단계: Bias / RMSE")

rmse_long <- conv_df %>%
  select(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean,
         a_est, b1_est, b2_est, b3_est, b4_est,
         a_true, b1_true, b2_true, b3_true, b4_true,
         err_a, err_b1, err_b2, err_b3, err_b4)
write.csv(rmse_long, out("rmse_bias_long.csv"), row.names = FALSE)

rmse_summary <- rmse_long %>%
  pivot_longer(starts_with("err_"), names_to = "param", values_to = "err",
               names_prefix = "err_") %>%
  group_by(cond_code, param, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
  summarise(n = sum(!is.na(err)),
            bias   = mean(err, na.rm = TRUE),
            rmse   = sqrt(mean(err^2, na.rm = TRUE)),
            sd_err = sd(err, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(cond_code, param)
write.csv(rmse_summary, out("rmse_bias_summary.csv"), row.names = FALSE)
lg("  저장 완료: rmse_bias_long.csv / rmse_bias_summary.csv")

b_params <- c("b1", "b2", "b3", "b4")
plot_b <- rmse_summary %>%
  filter(param %in% b_params) %>%
  add_factors() %>%
  mutate(param = factor(param, levels = b_params,
                        labels = c("b₁", "b₂", "b₃", "b₄")))

common_layers <- list(
  facet_grid(rows = vars(param), cols = vars(f_iv2, f_iv1),
             labeller = labeller(f_iv2 = function(x) paste0("Interval = ", x),
                                 f_iv1 = function(x) paste0("s4 = ", x),
                                 param = label_value)),
  labs(x = "b-parameter Mean (IV3)"),
  theme_bw(base_size = 14),
  theme(strip.text = element_text(size = 11), axis.text = element_text(size = 11))
)

p_bias <- ggplot(plot_b, aes(x = f_iv3, y = bias, group = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7, colour = "steelblue") +
  geom_point(size = 1.8, colour = "steelblue") +
  labs(y = "Bias") + common_layers
ggsave(out("rmse_bias_plot_b_bias.png"), p_bias, width = 16, height = 10, dpi = 150)

p_rmse <- ggplot(plot_b, aes(x = f_iv3, y = rmse, group = 1)) +
  geom_line(linewidth = 0.7, colour = "firebrick") +
  geom_point(size = 1.8, colour = "firebrick") +
  labs(y = "RMSE") + common_layers
ggsave(out("rmse_bias_plot_b_rmse.png"), p_rmse, width = 16, height = 10, dpi = 150)
lg("  저장 완료: rmse_bias_plot_b_bias.png / rmse_bias_plot_b_rmse.png")

# =============================================================================
# 5단계: ANOVA — b3 / b4 Bias·RMSE (Type III, IV1 * IV2 * IV3)
# =============================================================================
lg_section("5단계: ANOVA (b3/b4 Bias & RMSE)")

anova_data <- add_factors(rmse_long)
contr_list <- list(f_iv1 = contr.sum, f_iv2 = contr.sum, f_iv3 = contr.sum)

fit_type3_anova <- function(dv, label) {
  df <- anova_data %>% mutate(.dv = dv) %>% filter(!is.na(.dv))
  fit <- lm(.dv ~ f_iv1 * f_iv2 * f_iv3, data = df, contrasts = contr_list)
  tbl <- as.data.frame(car::Anova(fit, type = 3))
  tbl[["Mean Sq"]] <- tbl[["Sum Sq"]] / tbl[["Df"]]
  ss_tot <- sum(tbl[rownames(tbl) != "(Intercept)", "Sum Sq"], na.rm = TRUE)
  tbl$eta_sq <- tbl[["Sum Sq"]] / ss_tot
  lg(sprintf("  [%s] 완료 (n=%d)", label, nrow(df)))
  list(tbl = tbl, n = nrow(df))
}

anova_res <- list(
  "b3 — Bias" = fit_type3_anova(anova_data$err_b3,   "b3 Bias"),
  "b3 — RMSE" = fit_type3_anova(anova_data$err_b3^2, "b3 RMSE"),
  "b4 — Bias" = fit_type3_anova(anova_data$err_b4,   "b4 Bias"),
  "b4 — RMSE" = fit_type3_anova(anova_data$err_b4^2, "b4 RMSE")
)
dv_desc <- c("err_b3 = b3_est - b3_true (범주3-4 경계)", "err_b3^2",
             "err_b4 = b4_est - b4_true (범주4-5 경계)", "err_b4^2")

ttest_df <- bind_rows(lapply(c("b3", "b4"), function(p) {
  rmse_long %>%
    select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean,
           err = all_of(paste0("err_", p))) %>%
    filter(!is.na(err)) %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    filter(n() >= 2) %>%
    summarise(n = n(), mean_bias = mean(err),
              t_stat  = unname(t.test(err, mu = 0)$statistic),
              p_value = t.test(err, mu = 0)$p.value,
              .groups = "drop") %>%
    mutate(param = p, .before = 1)
}))
ttest_df$p_adj_fdr <- p.adjust(ttest_df$p_value, method = "BH")
ttest_df$sig <- ifelse(ttest_df$p_adj_fdr < 0.001, "***",
                ifelse(ttest_df$p_adj_fdr < 0.01,  "**",
                ifelse(ttest_df$p_adj_fdr < 0.05,  "*", "")))
write.csv(ttest_df, out("bias_ttest_by_cond.csv"), row.names = FALSE)

make_anova_df <- function(tbl) {
  data.frame(Df = tbl[, "Df"],
             SS = round(tbl[, "Sum Sq"], 4),
             MS = round(tbl[, "Mean Sq"], 6),
             F_value = round(tbl[, "F value"], 4),
             p_value = format.pval(tbl[, "Pr(>F)"], digits = 3, eps = 0.001),
             eta_sq  = round(tbl[, "eta_sq"], 4),
             row.names = rownames(tbl))
}

sink(out("bias_rmse_anova.txt"))
tryCatch({
  for (i in seq_along(anova_res)) {
    cat(strrep("=", 64), "\n")
    cat(sprintf("ANOVA [%s]  (IV4 = 정규분포 고정)\n종속변수: %s\n",
                names(anova_res)[i], dv_desc[i]))
    cat("독립변수: IV1 * IV2 * IV3  |  Type III SS (car::Anova, contr.sum)\n")
    cat(strrep("=", 64), "\n\n")
    cat(sprintf("수렴 성공 반복 수: %d\n\n", anova_res[[i]]$n))
    print(make_anova_df(anova_res[[i]]$tbl))
    cat("\n")
  }
  cat(strrep("=", 64), "\n")
  cat("조건별 단일표본 t검정 (mu=0) | FDR(BH) 보정\n")
  cat(strrep("=", 64), "\n\n")
  print(as.data.frame(ttest_df), row.names = FALSE, digits = 4)
}, error = function(e) {
  cat("\n[오류] 결과 저장 중 에러:", conditionMessage(e), "\n")
}, finally = sink())
lg("  저장 완료: bias_rmse_anova.txt / bias_ttest_by_cond.csv")

lg_section("완료")
lg(sprintf("  모든 결과 저장 위치: %s/", OUT_DIR))
flush_log()
