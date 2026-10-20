# =============================================================================
# 00_functions.R  (simulation1/R_graph)
# 그래프 스크립트 공통 함수: 패키지, 테마, 저장, 모형 확률, 집계
# 설정값은 00_config.R 에서 바꿈 (이 파일은 고치지 않아도 됨)
# =============================================================================

if (!file.exists("R_graph/00_config.R"))
  stop("작업 디렉토리를 simulation1/ 로 맞춘 뒤 실행하세요. 현재: ", getwd())
source("R_graph/00_config.R", encoding = "UTF-8")

# ── 패키지 ──────────────────────────────────────────────────────────────────
for (pkg in c("dplyr", "tidyr", "ggplot2", "scales", "ragg", "car")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.rstudio.com/")
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

for (d in c(DATA_DIR, FIG_DIR, TAB_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── 색 (색맹 친화 Okabe-Ito 계열) ────────────────────────────────────────────
COL_INTERVAL <- c("0.5" = "#E69F00", "1.0" = "#009E73", "1.5" = "#0072B2")
SHP_INTERVAL <- c("0.5" = 16, "1.0" = 17, "1.5" = 15)
LTY_INTERVAL <- c("0.5" = "solid", "1.0" = "dashed", "1.5" = "dotdash")
COL_PAIR     <- c("12" = "#56B4E9", "23" = "#E69F00", "34" = "#D55E00")
COL_CATEGORY <- c("1" = "#0072B2", "2" = "#56B4E9", "3" = "#009E73",
                  "4" = "#E69F00", "5" = "#D55E00")
COL_PARAM    <- c(b1 = "#56B4E9", b2 = "#009E73", b3 = "#0072B2", b4 = "#D55E00")
COL_TRUE     <- "black"
COL_GPCM     <- "#D55E00"
COL_BAR      <- "#4A7FB5"   # 전치 비율 막대
COL_BAR_NEG  <- "#C8553D"   # 차이 그래프의 음수 막대
COL_BELOW    <- "#1B9E3E"   # 경계모수보다 작은 θ (초록)
COL_ABOVE    <- "#1F5FD6"   # 경계모수보다 큰 θ (파랑)

# ── 테마 ────────────────────────────────────────────────────────────────────
theme_thesis <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size, base_family = FONT_FAMILY) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = "grey40"),
      strip.text        = element_text(size = rel(0.9), margin = margin(6, 6, 6, 6)),
      axis.title        = element_text(size = rel(1)),
      axis.text         = element_text(size = rel(0.85), colour = "black"),
      legend.position   = "bottom",
      legend.title      = element_text(size = rel(0.9)),
      legend.text       = element_text(size = rel(0.85)),
      legend.key.width  = unit(3, "lines"),
      plot.margin       = margin(10, 15, 10, 10)
    )
}

# 선·점 크기 (글씨 크기에 맞춤)
LW <- 1.3   # 선 굵기
PS <- 4     # 점 크기

# ── 저장 ────────────────────────────────────────────────────────────────────
save_fig <- function(p, name, width, height) {
  path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(path, p, width = width, height = height, units = "in", dpi = DPI,
         device = ragg::agg_png, bg = "white")
  if (isTRUE(SAVE_PDF))
    ggsave(sub("\\.png$", ".pdf", path), p, width = width, height = height,
           units = "in", device = grDevices::cairo_pdf, family = FONT_FAMILY)
  message("  그림 저장: ", path)
  invisible(path)
}

# 엑셀에서 한글이 깨지지 않도록 UTF-8 BOM 을 붙여 CSV 저장
save_table <- function(df, name) {
  path <- file.path(TAB_DIR, paste0(name, ".csv"))
  con  <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con))
  cat("\ufeff", file = con)
  op <- options(scipen = 999); on.exit(options(op), add = TRUE)
  write.csv(df, con, row.names = FALSE, na = "")
  message("  표 저장: ", path)
  invisible(path)
}

# ── 수준 라벨·factor ────────────────────────────────────────────────────────
fmt_level <- function(x) ifelse(abs(x) < 1e-9, "0", sprintf("%.1f", x))

SF_LEVELS       <- unname(SF_GAP_LABEL[as.character(sort(SIM1_IV1))])
INTERVAL_LEVELS <- fmt_level(sort(SIM1_IV2))
MEAN_LEVELS     <- fmt_level(sort(SIM1_IV3))

add_design_factors <- function(df) {
  df %>% mutate(
    f_sf   = factor(unname(SF_GAP_LABEL[as.character(sf4)]), levels = SF_LEVELS),
    f_int  = factor(fmt_level(interval), levels = INTERVAL_LEVELS),
    f_mean = factor(fmt_level(mean),     levels = MEAN_LEVELS)
  )
}

lab_sf       <- function(x) paste0(LAB$sf, " = ", x)
lab_interval <- function(x) paste0(LAB$interval, " = ", x)
lab_mean     <- function(x) paste0(LAB$mean, " = ", x)

# 조건 c(sf4, interval, mean) 와 일치하는 행 (실수 비교 허용 오차)
match_cond <- function(df, cond) {
  abs(df$sf4 - cond[["sf4"]]) < 1e-6 &
    abs(df$interval - cond[["interval"]]) < 1e-6 &
    abs(df$mean - cond[["mean"]]) < 1e-6
}
cond_text <- function(cond, sep = " · ") {
  paste(lab_sf(SF_GAP_LABEL[[as.character(cond[["sf4"]])]]),
        lab_interval(fmt_level(cond[["interval"]])),
        lab_mean(fmt_level(cond[["mean"]])), sep = sep)
}

# 시뮬레이션1 조건 코드 (정규분포)
sim1_code <- function(sf4, interval, mean) {
  paste0(match(sf4, SIM1_IV1), match(interval, SIM1_IV2),
         match(mean, SIM1_IV3), SIM1_NORMAL_DIGIT)
}

# ── 모형 확률 ───────────────────────────────────────────────────────────────
# 생성 모형 (채점함수 일반화 GPCM, simulation1/R/03_response.R 과 같은 식)
#   η_k = Σ_{j≤k} a·w_j·(θ − b_j),  w_j = s_j − s_{j−1}
#   반환: 행 = θ, 열 = 범주 1~5
#   = NRM 꼴로 쓰면 기울기 a·s_k, 절편 −a·Σ_{j≤k} w_j·b_j
crc_true <- function(theta, b, sf = 0:4, a = DISCRIM) {
  w <- diff(sf)
  crc_nrm(theta, a, c(0, cumsum(w)), -a * c(0, cumsum(w * b)))
}

# NRM (mirt nominal): η_k = a1·ak_k·θ + d_k,  k = 0..4
crc_nrm <- function(theta, a1, ak, d) {
  eta <- outer(theta, a1 * ak) + matrix(d, length(theta), length(d), byrow = TRUE)
  eta <- eta - apply(eta, 1, max)
  p   <- exp(eta)
  p / rowSums(p)
}

crc_long <- function(theta, P, ...) {
  data.frame(theta = rep(theta, ncol(P)),
             category = factor(rep(seq_len(ncol(P)), each = length(theta))),
             prob = as.vector(P), ...)
}

# 범주 번호를 각 곡선의 가장 높은 곳에 적기 위한 위치
crc_peak_labels <- function(df, group_vars = character(0)) {
  df %>%
    group_by(across(all_of(c(group_vars, "category")))) %>%
    slice_max(prob, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(theta = pmin(pmax(theta, min(df$theta) + 0.25), max(df$theta) - 0.25))
}

# ── 자료 불러오기 (01_collect_data.R 결과) ──────────────────────────────────
load_sim1 <- function() {
  p <- file.path(DATA_DIR, "sim1_gpcm_item1.rds")
  if (!file.exists(p)) stop("먼저 R_graph/01_collect_data.R 을 실행하세요: ", p)
  readRDS(p)
}
load_sim2 <- function() {
  p <- file.path(DATA_DIR, "sim2_nrm_item1.rds")
  if (!file.exists(p)) stop("먼저 R_graph/01_collect_data.R 을 실행하세요: ", p)
  readRDS(p)
}

# 수렴하고 b1~b4 가 모두 유한한 반복 + 오차·전치 표시
prepare_reps <- function(df) {
  df %>%
    filter(converged %in% TRUE,
           is.finite(b1), is.finite(b2), is.finite(b3), is.finite(b4)) %>%
    mutate(
      err_b1 = b1 - b1_true, err_b2 = b2 - b2_true,
      err_b3 = b3 - b3_true, err_b4 = b4 - b4_true,
      t12 = b1 >= b2, t23 = b2 >= b3, t34 = b3 >= b4,
      transposed = t12 | t23 | t34
    ) %>%
    add_design_factors()
}

# NRM 이상치·채점함수 이상 표시
flag_nrm <- function(df) {
  dak <- cbind(df$ak1 - df$ak0, df$ak2 - df$ak1, df$ak3 - df$ak2, df$ak4 - df$ak3)
  err <- abs(cbind(df$err_b1, df$err_b2, df$err_b3, df$err_b4))
  df$extreme  <- apply(err >= ERR_CUT, 1, any)
  df$sf_issue <- apply(dak <= 0 | abs(dak) < DAK_CUT, 1, any)
  df$sf_rev4  <- dak[, 4] <= 0
  df
}

# 조건별 전치 비율 (분모: 수렴 + 경계모수 계산된 반복)
summarise_transposition <- function(reps) {
  reps %>%
    group_by(f_sf, f_int, f_mean) %>%
    summarise(n = n(),
              prop_any = mean(transposed),
              prop_12  = mean(t12), prop_23 = mean(t23), prop_34 = mean(t34),
              .groups = "drop")
}

# 조건 × 경계모수별 편향·RMSE
#   exclude: 오차 계산에서 뺄 값 표시 (반복 × 경계모수 논리 행렬, NULL = 모두 사용)
summarise_error <- function(reps, exclude = NULL) {
  long <- reps %>%
    select(f_sf, f_int, f_mean, err_b1, err_b2, err_b3, err_b4) %>%
    pivot_longer(starts_with("err_"), names_to = "param",
                 names_prefix = "err_", values_to = "err")
  if (!is.null(exclude)) long$err[as.vector(t(exclude))] <- NA
  long %>%
    group_by(f_sf, f_int, f_mean, param) %>%
    summarise(n = sum(!is.na(err)),
              bias = mean(err, na.rm = TRUE),
              rmse = sqrt(mean(err^2, na.rm = TRUE)),
              .groups = "drop")
}

# NRM 이상치 제거 표시 (OUTLIER_UNIT 에 따라 반복 전체 또는 해당 값만)
nrm_outlier_matrix <- function(reps) {
  err <- abs(cbind(reps$err_b1, reps$err_b2, reps$err_b3, reps$err_b4)) >= ERR_CUT
  if (OUTLIER_UNIT == "rep") err[] <- apply(err, 1, any)
  err
}

# ── 효과 이름 (표용) ────────────────────────────────────────────────────────
EFFECT_LABEL <- c(f_sf = LAB$sf, f_int = LAB$interval, f_mean = LAB$mean)
effect_name <- function(term) {
  vapply(strsplit(term, ":"), function(parts) {
    lab <- EFFECT_LABEL[parts]
    lab[is.na(lab)] <- parts[is.na(lab)]
    paste(lab, collapse = " × ")
  }, character(1))
}
fmt_p <- function(p) ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p)))

# ── 공통 y 범위 ─────────────────────────────────────────────────────────────
# 여러 자료의 값을 모두 담는 보기 좋은 범위 (비교할 그림끼리 같은 값을 넘겨 축을 맞춤)
common_ylim <- function(values, include_zero = TRUE) {
  v <- values[is.finite(values)]
  range(pretty(c(v, if (include_zero) 0)))
}
nice_upper <- function(x, step = 0.1) ceiling(max(x, na.rm = TRUE) / step) * step

# 9개 패널: 행 = 채점함수 간격, 열 = 경계모수 간격
facet9 <- function() {
  facet_grid(rows = vars(f_sf), cols = vars(f_int),
             labeller = labeller(f_sf = lab_sf, f_int = lab_interval))
}
pct_labels <- scales::label_percent(accuracy = 1)

# ── 공통 그래프 1: 조건별 선 그래프 (편향·RMSE) ─────────────────────────────
#   x = 문항 난이도, 색·모양·선 = 경계모수 간격, 열 = 채점함수 간격, (행 = rows 변수)
#   ylim: 모든 패널이 같은 y축 (비교할 그림끼리 같은 ylim 을 넘김)
plot_lines <- function(df, y, ylab, ylim, rows = NULL, hline = NULL,
                       y_labels = waiver(), row_labeller = label_parsed) {
  p <- ggplot(df, aes(f_mean, .data[[y]], colour = f_int, shape = f_int,
                      linetype = f_int, group = f_int))
  if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, colour = "grey45", linewidth = 0.7)
  p <- p +
    geom_line(linewidth = LW) +
    geom_point(size = PS) +
    scale_colour_manual(values = COL_INTERVAL, name = LAB$interval, drop = FALSE) +
    scale_shape_manual(values = SHP_INTERVAL, name = LAB$interval, drop = FALSE) +
    scale_linetype_manual(values = LTY_INTERVAL, name = LAB$interval, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = y_labels, limits = ylim) +
    labs(x = LAB$mean, y = ylab) +
    theme_thesis()
  if (is.null(rows)) {
    p + facet_grid(cols = vars(f_sf), labeller = labeller(f_sf = lab_sf))
  } else {
    p + facet_grid(rows = vars(.data[[rows]]), cols = vars(f_sf),
                   labeller = labeller(.rows = row_labeller, f_sf = lab_sf)) +
      theme(strip.text.y = element_text(angle = 0, size = rel(1.15)))
  }
}

# ── 공통 그래프 2: 9개 패널 막대 (전치 비율) ────────────────────────────────
#   행 = 채점함수 간격, 열 = 경계모수 간격, x = 문항 난이도
#   fill_by_sign = TRUE 면 양수·음수 막대 색을 다르게 (차이 그래프)
plot_bars <- function(df, y, ylab, ylim, y_labels = pct_labels, fill_by_sign = FALSE) {
  p <- ggplot(df, aes(f_mean, .data[[y]]))
  if (fill_by_sign) {
    p <- p + geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.6) +
      geom_col(aes(fill = .data[[y]] >= 0), width = 0.7, colour = "grey30", linewidth = 0.3) +
      scale_fill_manual(values = c(`TRUE` = COL_BAR, `FALSE` = COL_BAR_NEG), guide = "none")
  } else {
    p <- p + geom_col(fill = COL_BAR, width = 0.7, colour = "grey30", linewidth = 0.3)
  }
  p + facet9() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = y_labels, limits = ylim) +
    labs(x = LAB$mean, y = ylab) +
    theme_thesis()
}

# ── 공통 그래프 3: 9개 패널 누적 막대 ────────────────────────────────────────
#   plot_bars 와 같은 배치, 막대 = fill 변수 누적
plot_stack <- function(df, fill, values, labels, legend_title, ylim,
                       total = NULL, total_label = NULL) {
  p <- ggplot(df, aes(f_mean, prop)) +
    geom_col(aes(fill = .data[[fill]]), position = position_stack(reverse = TRUE),
             width = 0.7, colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = values, labels = labels, name = legend_title, drop = FALSE)
  if (!is.null(total)) {
    p <- p + geom_point(data = total, aes(y = prop_any, shape = "total"), size = PS + 1,
                        colour = "black") +
      scale_shape_manual(values = c(total = 18), labels = total_label, name = NULL)
  }
  p + facet9() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = pct_labels, limits = ylim) +
    labs(x = LAB$mean, y = LAB$trans) +
    theme_thesis() +
    theme(legend.box = "vertical", legend.key.width = unit(1.5, "lines"),
          legend.key.height = unit(1.6, "lines"))
}

# ── 로지스틱 회귀 우도비 검정 표 (전치 여부 ~ 채점함수 × 경계모수 간격 × 난이도) ──
transposition_lrt <- function(reps) {
  d   <- reps %>% mutate(y = as.numeric(transposed))
  fit <- suppressWarnings(glm(y ~ f_sf * f_int * f_mean, data = d, family = binomial("logit")))
  lrt <- suppressWarnings(drop1(fit, scope = ~., test = "Chisq"))
  col <- intersect(c("LRT", "Chisq"), colnames(lrt))[1]
  tm  <- rownames(lrt)[rownames(lrt) != "<none>"]
  data.frame(`효과` = effect_name(tm),
             `자유도` = lrt[tm, "Df"],
             `우도비 카이제곱` = round(lrt[tm, col], 2),
             p = fmt_p(lrt[tm, grep("^Pr", colnames(lrt))[1]]),
             `편 의사 R제곱` = round(lrt[tm, col] / fit$null.deviance, 4),
             check.names = FALSE)
}

# ── Type III 분산분석 표 (오차 ~ 채점함수 × 경계모수 간격 × 난이도) ──────────
anova_type3 <- function(reps, dv) {
  d   <- reps %>% mutate(.dv = .data[[dv]]) %>% filter(is.finite(.dv))
  fit <- lm(.dv ~ f_sf * f_int * f_mean, data = d,
            contrasts = list(f_sf = contr.sum, f_int = contr.sum, f_mean = contr.sum))
  tb  <- as.data.frame(car::Anova(fit, type = 3))
  tb  <- tb[rownames(tb) != "(Intercept)", ]
  ss_res <- tb["Residuals", "Sum Sq"]
  ss_tot <- sum(tb[["Sum Sq"]])
  data.frame(`효과` = ifelse(rownames(tb) == "Residuals", "오차", effect_name(rownames(tb))),
             `제곱합` = round(tb[["Sum Sq"]], 3),
             `자유도` = tb[["Df"]],
             `평균제곱` = round(tb[["Sum Sq"]] / tb[["Df"]], 4),
             F = round(tb[["F value"]], 2),
             p = ifelse(is.na(tb[["Pr(>F)"]]), "", fmt_p(tb[["Pr(>F)"]])),
             `에타제곱` = ifelse(rownames(tb) == "Residuals", NA, round(tb[["Sum Sq"]] / ss_tot, 4)),
             `편 에타제곱` = ifelse(rownames(tb) == "Residuals", NA,
                                round(tb[["Sum Sq"]] / (tb[["Sum Sq"]] + ss_res), 4)),
             check.names = FALSE)
}

# 전치 비율 막대 그림들(그림 7·8·12·15, 예비 NRM 쌍별)이 함께 쓰는 y 범위
#   쌍별 누적 막대는 한 반복에서 여러 쌍이 전치되면 전체 비율보다 높아질 수 있어 그 값까지 포함
trans_ylim <- function(sums) {
  c(0, nice_upper(unlist(lapply(sums, function(s)
    pmax(s$prop_any, s$prop_12 + s$prop_23 + s$prop_34)))))
}
