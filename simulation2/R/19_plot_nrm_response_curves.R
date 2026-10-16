# =============================================================================
# 19_plot_nrm_response_curves.R  (simulation2, NRM 추정 결과)
# 문항 1의 생성 모형 범주반응곡선과, 추정된 NRM 모수의 평균(또는 중앙값)으로 만든
# 범주반응곡선을 겹쳐 그림
#
# ── 범주반응곡선 (category response curve) ───────────────────────────────────
#   P(X=k | θ),  k = 0..4  (그래프 표기는 범주 1~5)
#     생성 모형(GPCM, 채점함수 0,1,2,s4,4): η_k = Σ_{j≤k} a · w_j · (θ − b_j),  w_j = sf_j − sf_{j−1}
#     NRM 추정                            : η_k = α · s_k · θ + c_k
#     P(X=k | θ) = exp(η_k) / Σ_m exp(η_m)
#   생성 모형을 NRM 모수로 쓰면 α = 1, s = (0, 1, 2, s4, 4), c_k = −Σ_{j≤k} w_j b_j
#
# ── 추정 곡선의 모수 ─────────────────────────────────────────────────────────
#   조건마다 반복들의 α, s1~s3, c1~c4 를 SUMMARY_STAT("mean" 또는 "median")으로 요약해
#   곡선 1개를 만듦 (반복마다 곡선을 그려 평균한 것이 아님)
#   반복 집합은 16·18번과 같은 3가지: all / excl_extreme / excl_extreme_sfrev
#   ※ 평균은 극단 반복에 크게 흔들리므로 all 집합은 중앙값이나 제외 집합과 함께 볼 것
#
# 입력: output/nrm/analysis/nrm_test_rep.csv  (12_analysis_nrm_test.R 출력, 문항 1 모수만 있음)
# 출력: output/nrm/response_curves/<반복 집합>/
#   curve_params_<stat>.csv          — 조건별 요약 모수(α, s0~s4, c0~c4)와 그로 계산한 경계모수 b1~b4
#   b_comparison_<stat>.csv          — 조건 × b1~b4: 참값, 요약 모수로 계산한 b(곡선의 경계),
#                                      b̂ 의 평균·중앙값, 각각의 편향, Δs_k ≤ 0 인 반복 비율
#     ※ Bias(18·21번 등) = mean(b̂) − b 이고, 곡선의 경계는 mean(α, s, c) 로 계산한 b 라서 서로 다름
#       b_k = (c_{k−1} − c_k)/(α Δs_k) 가 비율이므로, Δs_k 가 0 근처이거나 음수(전치)인 반복에서
#       b̂_k 가 크게 튀거나 부호가 뒤집혀 b̂ 평균만 끌려감 → 곡선은 참값과 비슷해도 Bias 는 클 수 있음
#   crc_<stat>_mean<IV3>.png         — 범주반응곡선 (IV3 수준마다 1장, 행 IV1 × 열 IV2)
#   생성 모형 = 점선, NRM 추정 = 실선, 색 = 범주(1~5), 회색 세로선 = 생성 모형의 경계모수 b1~b4
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
SUMMARY_STAT <- c("mean", "median")  # 둘 다 그리려면 c("mean", "median")
THETA_RANGE  <- c(-4, 5)
ERR_CUT      <- 3     # 극단 오차: |b̂ − b| ≥ ERR_CUT   (15·16·18번과 같은 값)
DAK_CUT      <- 0.1   # 분모 붕괴: |Δs_k| < DAK_CUT      (15·16·18번과 같은 값)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

in_path <- "output/nrm/analysis/nrm_test_rep.csv"
if (!file.exists(in_path)) stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", in_path)
BASE_DIR <- "output/nrm/response_curves"
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)
THETA    <- seq(THETA_RANGE[1], THETA_RANGE[2], length.out = 301)

# =============================================================================
# 1. 데이터 (16·18번과 같은 제외 기준)
# =============================================================================
rep_df <- read.csv(in_path, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE,
         !is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4))

bt <- sapply(B_OFFSET, function(o) rep_df$iv3_b_mean + rep_df$iv2_b_interval * o)
rep_df <- rep_df %>%
  mutate(is_extreme = pmax(abs(b1 - bt[, 1]), abs(b2 - bt[, 2]),
                           abs(b3 - bt[, 3]), abs(b4 - bt[, 4])) >= ERR_CUT,
         is_sfrev   = {
           dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
           rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
         })

# =============================================================================
# 2. 곡선 계산
# =============================================================================
# 생성 모형의 NRM 표현: α = 1, s = (0,1,2,s4,4), c_k = −Σ_{j≤k} w_j b_j
true_nrm <- function(sf4, b_int, b_mean) {
  s <- c(0, 1, 2, sf4, 4)
  b <- b_mean + b_int * B_OFFSET
  list(alpha = 1, s = s, c = c(0, -cumsum(diff(s) * b)))
}
# 범주반응곡선: 행 θ, 열 범주 0~4
crc <- function(alpha, s, c) {
  z <- outer(THETA, alpha * s) + matrix(c, length(THETA), 5, byrow = TRUE)
  z <- z - apply(z, 1, max); p <- exp(z); p / rowSums(p)
}
curve_long <- function(mat, names, model) {
  colnames(mat) <- names
  data.frame(theta = THETA, mat, check.names = FALSE) %>%
    pivot_longer(-theta, names_to = "curve", values_to = "p") %>%
    mutate(model = model)
}

# =============================================================================
# 3. 그래프
# =============================================================================
make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
MODEL_LEVELS <- c("Generating model", "NRM estimate")
CRC_COLOURS  <- c(`Cat 1` = "#4C9BE8", `Cat 2` = "#59A14F", `Cat 3` = "#F2A541",
                  `Cat 4` = "#D1495B", `Cat 5` = "#6A3D9A")

curve_plot <- function(df, vlines, colours, ylab, title, path) {
  p <- ggplot(df, aes(x = theta, y = p, colour = curve, linetype = model,
                      group = interaction(curve, model))) +
    geom_vline(data = vlines, aes(xintercept = b), colour = "grey75", linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_colour_manual(values = colours, name = NULL) +
    scale_linetype_manual(values = c("Generating model" = "dashed", "NRM estimate" = "solid"),
                          name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(title = title, x = expression(theta), y = ylab) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 7), axis.text = element_text(size = 7))
  ggsave(path, p, width = 24, height = 18, dpi = 150)
}

run_set <- function(keep, label, desc) {
  out_dir <- file.path(BASE_DIR, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)
  kept <- rep_df[keep, ]

  for (stat in SUMMARY_STAT) {
    fun <- switch(stat, mean = mean, median = median, stop("SUMMARY_STAT: mean 또는 median"))

    # 조건별 요약 모수 (s0 = 0, s4 = 4, c0 = 0 은 식별 제약으로 고정)
    par_df <- kept %>%
      group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
      summarise(n = n(),
                alpha = fun(a1),
                s0 = fun(ak0), s1 = fun(ak1), s2 = fun(ak2), s3 = fun(ak3), s4 = fun(ak4),
                c0 = fun(d0), c1 = fun(d1), c2 = fun(d2), c3 = fun(d3), c4 = fun(d4),
                .groups = "drop") %>%
      mutate(b1_curve = (c0 - c1) / (alpha * (s1 - s0)),
             b2_curve = (c1 - c2) / (alpha * (s2 - s1)),
             b3_curve = (c2 - c3) / (alpha * (s3 - s2)),
             b4_curve = (c3 - c4) / (alpha * (s4 - s3)))
    write.csv(par_df, out(sprintf("curve_params_%s.csv", stat)), row.names = FALSE)

    # 요약 모수로 계산한 경계모수(곡선의 경계) vs 반복별 경계모수 추정값의 평균·중앙값
    # b_k = (c_{k-1} − c_k) / (α Δs_k) 는 비선형(비율)이라 "모수 평균으로 만든 b" ≠ "b̂ 의 평균"
    # 특히 Δs_k ≤ 0 (채점함수 전치) 인 반복은 b̂_k 부호가 뒤집혀 b̂ 평균만 크게 끌어당김
    b_cmp <- bind_rows(lapply(1:4, function(k) {
      dk <- kept[[paste0("ak", k)]] - kept[[paste0("ak", k - 1)]]
      data.frame(cond_code = kept$cond_code, param = paste0("b", k),
                 b_hat = kept[[paste0("b", k)]], dak = dk, stringsAsFactors = FALSE)
    })) %>%
      filter(!is.nan(b_hat)) %>%
      group_by(cond_code, param) %>%
      summarise(mean_b_hat   = mean(b_hat[is.finite(b_hat)]),
                median_b_hat = median(b_hat),
                share_ds_le0 = mean(dak <= 0),
                .groups = "drop") %>%
      left_join(par_df %>%
                  select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, starts_with("b") & ends_with("_curve")) %>%
                  pivot_longer(ends_with("_curve"), names_to = "param", values_to = "b_curve") %>%
                  mutate(param = sub("_curve$", "", param)),
                by = c("cond_code", "param")) %>%
      mutate(b_true = iv3_b_mean + iv2_b_interval * B_OFFSET[as.integer(sub("b", "", param))],
             bias_mean_b_hat   = mean_b_hat - b_true,
             bias_median_b_hat = median_b_hat - b_true,
             bias_b_curve      = b_curve - b_true) %>%
      select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param, b_true, b_curve,
             mean_b_hat, median_b_hat, bias_b_curve, bias_mean_b_hat, bias_median_b_hat, share_ds_le0) %>%
      arrange(cond_code, param)
    write.csv(b_cmp, out(sprintf("b_comparison_%s.csv", stat)), row.names = FALSE)

    curves <- bind_rows(lapply(seq_len(nrow(par_df)), function(i) {
      r  <- par_df[i, ]
      tp <- true_nrm(r$iv1_sf4, r$iv2_b_interval, r$iv3_b_mean)
      es <- c(r$s0, r$s1, r$s2, r$s3, r$s4)
      ec <- c(r$c0, r$c1, r$c2, r$c3, r$c4)
      key <- r[, c("cond_code", "iv1_sf4", "iv2_b_interval", "iv3_b_mean")]
      bind_rows(
        curve_long(crc(tp$alpha, tp$s, tp$c), paste("Cat", 1:5), MODEL_LEVELS[1]),
        curve_long(crc(r$alpha, es, ec),      paste("Cat", 1:5), MODEL_LEVELS[2])
      ) %>% cbind(key, .)
    })) %>%
      mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
             iv2 = make_num_factor(iv2_b_interval, "%.2f"),   # 좁음 → 넓음
             model = factor(model, levels = MODEL_LEVELS),
             curve = factor(curve, levels = names(CRC_COLOURS)))
    vlines <- par_df %>%
      select(iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
      tidyr::crossing(offset = B_OFFSET) %>%
      mutate(b   = iv3_b_mean + iv2_b_interval * offset,
             iv1 = make_num_factor(iv1_sf4,        "%.2f"),
             iv2 = make_num_factor(iv2_b_interval, "%.2f"))

    for (m in sort(unique(curves$iv3_b_mean))) {
      sub <- curves %>% filter(iv3_b_mean == m)
      tag <- sprintf("%s_mean%.1f", stat, m)
      curve_plot(sub, vlines %>% filter(iv3_b_mean == m), CRC_COLOURS, "P(X = k)",
                 sprintf("Item 1 category response curves, b-mean = %.1f (NRM parameters: %s) - %s", m, stat, desc),
                 out(sprintf("crc_%s.png", tag)))
    }
  }
  message(sprintf("[%s] 반복 %d개 → %s/", label, nrow(kept), out_dir))
}

run_set(rep(TRUE, nrow(rep_df)), "all", "all converged replications")
run_set(!rep_df$is_extreme, "excl_extreme",
        sprintf("excluding extreme errors (|error| >= %g)", ERR_CUT))
run_set(!rep_df$is_extreme & !rep_df$is_sfrev, "excl_extreme_sfrev",
        sprintf("excluding extreme errors and SF reversals / near-zero ds (|ds| < %g)", DAK_CUT))
