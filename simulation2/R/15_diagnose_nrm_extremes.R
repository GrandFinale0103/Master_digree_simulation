# =============================================================================
# 15_diagnose_nrm_extremes.R  (simulation2, NRM 추정 결과)
# NRM 경계모수 추정의 극단 오차 진단
#
# 판정 기준
#   극단 오차   : |b̂_k − b_k| ≥ ERR_CUT (3)                       — 연구자 설정 기준
#   부적절한 해 : Δak_k = ak_k − ak_{k−1} ≤ 0 (채점함수 전치)
#                 또는 |Δak_k| < DAK_CUT (0.1)                     — 경계모수 식의 분모 붕괴
#   (NRM 극단값 고정 기준을 제시한 문헌은 찾지 못함 → 다분 IRT 시뮬레이션 문헌의
#    "부적절한 해 판정 + 제외 수 보고 + 중앙값 기반 지표 병기" 원칙을 따름)
#
# 원인 표시 (극단 오차 경계모수마다, 중복 가능)
#   c1_near_zero : |Δak_k| < DAK_CUT
#   c2_sf_rev    : Δak_k ≤ 0
#   c3_sparse    : 경계모수 양옆 범주 중 하나의 응답 수 ≤ SPARSE_CUT (유령 응답자 제외)
#   c4_ghost     : 유령 응답자가 양옆 범주 중 하나를 채우려고 추가됨
#   c5_se_fail   : b4−b3 차이의 SE 계산 실패
#
# 진단 두 가지를 따로 실행
#   all_b : b1~b4 전체   → output/nrm/diagnostics/all_b/
#   b34   : b3, b4 만    → output/nrm/diagnostics/b34/
#
# 입력
#   output/nrm/analysis/nrm_test_rep.csv          (12_analysis_nrm_test.R)
#   output/analysis/resp_dist_rep_item1.csv       (07_analysis_response_dist.R, 유령 제외 응답 수)
#   output/responses/*_response.csv               (유령 응답자 정보: 극단 반복 + 정상 반복 표본)
#
# 출력 (진단 폴더마다)
#   extreme_rep.csv         — 반복 × 경계모수별 오차·판정·응답 수·분자/분모·원인
#   extreme_cond.csv        — 조건 × 경계모수별 극단·부적절 수, Bias/RMSE(전체·제외), 중앙값 편향, MdAE
#   extreme_cond_rep.csv    — 조건별 극단 반복 비율
#   extreme_cause_table.csv — 원인 조합별 빈도
#   extreme_summary.txt     — 기준, 제외 수, 극단×부적절 교차표, 원인 분포, 유령 비교
#   extreme_rate_plot.png   — 조건별 극단 반복 비율 (기존 전치 그래프 배치)
#   extreme_rate_by_param_bar.png — 조건별·경계모수별 극단 오차 비율 막대 (같은 배치)
#   extreme_resp_dist.png   — 극단/정상 반복의 범주 응답 비율 상자그림
#   extreme_vs_sparse.png   — 양옆 범주 최소 응답 수 vs |오차|
#   case_<조건>_rep<번호>_b<k>.png — 대표 사례: 참/추정 범주 반응 곡선
#   extreme_glm.txt         — 극단 반복 여부 ~ IV1*IV2*IV3 로지스틱 GLM (LRT)
#
# all_b 폴더 추가 출력
#   transposition_decomposition.png / .csv — 전치(인접 쌍 중 하나라도)를 원인별(둘 다/극단 오차만/
#                                            채점함수 이상만/둘 다 아님) 누적 막대로 분해 (분모: 수렴 반복)
#                                            극단 오차 = b1~b4 중 하나라도 |오차| ≥ ERR_CUT
#                                            채점함수 이상 = Δak1~Δak4 중 하나라도 Δak ≤ 0 또는 |Δak| < DAK_CUT
#   transposition_decomposition_share.png  — 전치된 반복 안에서의 원인 구성비 (작은 범주 확인용)
#   transposition_decomposition_counts.txt — 범주별 반복 수 합계
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
ERR_CUT          <- 3      # 극단 오차 기준 |오차| ≥ ERR_CUT
DAK_CUT          <- 0.1    # 부적절한 해: |Δak| < DAK_CUT
SPARSE_CUT       <- 5      # 희소 범주: 응답 수 ≤ SPARSE_CUT
EXCLUDE_IMPROPER <- TRUE   # 제외 지표 계산 시 부적절한 해도 제외할지 (FALSE면 극단 오차만 제외)
N_CASES          <- 6      # 대표 사례 수
NORMAL_SAMPLE    <- 2000   # 유령 응답자 비교용 정상 반복 표본 수
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

REP_PATH  <- "output/nrm/analysis/nrm_test_rep.csv"
DIST_PATH <- "output/analysis/resp_dist_rep_item1.csv"
RESP_DIR  <- "output/responses"
DIAG_DIR  <- "output/nrm/diagnostics"
if (!file.exists(REP_PATH))  stop("먼저 R/12_analysis_nrm_test.R 을 실행하세요: ", REP_PATH)
if (!file.exists(DIST_PATH)) stop("먼저 R/07_analysis_response_dist.R 을 실행하세요: ", DIST_PATH)

make_num_factor <- function(x, fmt, decreasing = FALSE) {
  vals <- sort(unique(na.omit(x)), decreasing = decreasing)
  factor(sprintf(fmt, x), levels = sprintf(fmt, vals))
}
add_factors <- function(df) {
  df %>% mutate(iv1 = make_num_factor(iv1_sf4,        "%.2f"),
                iv2 = make_num_factor(iv2_b_interval, "%.2f"),   # 좁은 → 넓은
                iv3 = make_num_factor(iv3_b_mean,     "%.1f"))
}

# =============================================================================
# 1. 데이터 결합
# =============================================================================
message("데이터 로딩 중...")
rep_df <- read.csv(REP_PATH, stringsAsFactors = FALSE,
                   colClasses = c(cond_code = "character")) %>%
  filter(as.logical(converged) %in% TRUE) %>%
  mutate(se_ok = as.logical(se_ok))

dist_df <- read.csv(DIST_PATH, stringsAsFactors = FALSE,
                    colClasses = c(cond_code = "character")) %>%
  select(cond_code, rep_id, n_persons, n0, n1, n2, n3, n4, prop0, prop1, prop2, prop3, prop4)

rep_df <- rep_df %>% left_join(dist_df, by = c("cond_code", "rep_id"))
if (any(is.na(rep_df$n0)))
  warning(sprintf("응답 분포가 없는 반복 %d개 (07번 결과와 반복이 다름)", sum(is.na(rep_df$n0))))

# ── 반복 × 경계모수 long 형식 ────────────────────────────────────────────────
long_df <- bind_rows(lapply(1:4, function(k) {
  # k 를 쓰는 계산은 transmute 밖에서 (안에서 k 열을 만들면 이후 k 가 열을 가리킴)
  sf       <- cbind(0, 1, 2, rep_df$iv1_sf4, 4)           # 참 채점함수 = 참 ak (a = 1)
  dak_vec  <- rep_df[[paste0("ak", k)]] - rep_df[[paste0("ak", k - 1)]]
  dakt_vec <- sf[, k + 1] - sf[, k]
  bt_vec   <- rep_df$iv3_b_mean + rep_df$iv2_b_interval * c(-1.5, -0.5, 0.5, 1.5)[k]
  bh_vec   <- rep_df[[paste0("b", k)]]
  num_vec  <- rep_df[[paste0("d", k - 1)]] - rep_df[[paste0("d", k)]]
  nlo_vec  <- rep_df[[paste0("n", k - 1)]]                # 범주 k   (1부터 센 번호)
  nhi_vec  <- rep_df[[paste0("n", k)]]                    # 범주 k+1
  rep_df %>% transmute(
    cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean, se_ok, a1,
    param     = paste0("b", k),
    b_hat     = bh_vec,
    b_true    = bt_vec,
    err       = bh_vec - bt_vec,
    dak       = dak_vec,
    dak_true  = dakt_vec,
    numer     = num_vec,
    denom     = a1 * dak_vec,
    n_lo      = nlo_vec,
    n_hi      = nhi_vec,
    n_min_adj = pmin(nlo_vec, nhi_vec),
    k         = k
  )
})) %>%
  mutate(
    extreme      = abs(err) >= ERR_CUT,
    improper     = dak <= 0 | abs(dak) < DAK_CUT,
    c1_near_zero = abs(dak) < DAK_CUT,
    c2_sf_rev    = dak <= 0,
    c3_sparse    = n_min_adj <= SPARSE_CUT,
    c5_se_fail   = !(se_ok %in% TRUE)
  )

# ── 유령 응답자 정보 (극단 반복 전부 + 정상 반복 표본) ──────────────────────
resp_files <- list.files(RESP_DIR, pattern = "_response\\.csv$", full.names = TRUE)
resp_index <- data.frame(
  path      = resp_files,
  cond_code = sub(".*_cond([0-9]+)_.*", "\\1", basename(resp_files)),
  rep_id    = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", basename(resp_files))),
  stringsAsFactors = FALSE
) %>% distinct(cond_code, rep_id, .keep_all = TRUE)

ext_keys <- long_df %>% filter(extreme) %>% distinct(cond_code, rep_id)
set.seed(20261011)
norm_all <- long_df %>%
  group_by(cond_code, rep_id) %>% filter(!any(extreme, na.rm = TRUE)) %>% ungroup() %>%
  distinct(cond_code, rep_id)
norm_keys <- norm_all %>% slice_sample(n = min(NORMAL_SAMPLE, nrow(norm_all)))
ghost_keys <- bind_rows(ext_keys, norm_keys) %>% inner_join(resp_index, by = c("cond_code", "rep_id"))

message(sprintf("유령 응답자 정보 읽는 중: 응답 파일 %d개", nrow(ghost_keys)))
ghost_info <- bind_rows(lapply(seq_len(nrow(ghost_keys)), function(i) {
  d <- read.csv(ghost_keys$path[i])
  g <- is.na(d$theta)
  data.frame(cond_code = ghost_keys$cond_code[i], rep_id = ghost_keys$rep_id[i],
             n_ghost = sum(g),
             ghost_cats = paste(sort(unique(d$item1[g])), collapse = ";"),  # 0부터 센 범주
             stringsAsFactors = FALSE)
}))
if (nrow(ghost_info) == 0)
  ghost_info <- data.frame(cond_code = character(0), rep_id = integer(0),
                           n_ghost = integer(0), ghost_cats = character(0))

long_df <- long_df %>%
  left_join(ghost_info, by = c("cond_code", "rep_id")) %>%
  mutate(c4_ghost = !is.na(ghost_cats) & mapply(function(gc, kk) {
    if (is.na(gc) || gc == "") return(FALSE)
    any(as.integer(strsplit(gc, ";")[[1]]) %in% c(kk - 1, kk))   # 양옆 범주(0부터): k−1, k
  }, ghost_cats, k))

# =============================================================================
# 2. 진단 함수 (대상 경계모수만 바꿔 두 번 실행)
# =============================================================================
nrm_crc <- function(theta, a1, ak, d) {
  z <- outer(theta, a1 * ak) + matrix(d, length(theta), 5, byrow = TRUE)
  z <- z - apply(z, 1, max); p <- exp(z); p / rowSums(p)
}
gpcm_crc <- function(theta, b, sf) {
  w   <- diff(sf)
  eta <- matrix(0, length(theta), 5)
  for (k in 2:5) eta[, k] <- eta[, k - 1] + w[k - 1] * (theta - b[k - 1])
  eta <- eta - apply(eta, 1, max); p <- exp(eta); p / rowSums(p)
}

run_diag <- function(params, label) {
  out_dir <- file.path(DIAG_DIR, label)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(f) file.path(out_dir, f)
  message(sprintf("\n=== 진단 [%s]: %s ===", label, paste(params, collapse = ", ")))

  ld <- long_df %>% filter(param %in% params)
  ld <- ld %>% mutate(excluded = extreme | (EXCLUDE_IMPROPER & improper))

  write.csv(ld %>% select(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean, param,
                          b_hat, b_true, err, extreme, improper, excluded,
                          dak, dak_true, numer, denom, n_lo, n_hi, n_min_adj,
                          n_ghost, ghost_cats, se_ok,
                          c1_near_zero, c2_sf_rev, c3_sparse, c4_ghost, c5_se_fail),
            out("extreme_rep.csv"), row.names = FALSE)

  # ── 조건 × 경계모수 지표 ──
  cond_df <- ld %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, param) %>%
    summarise(
      n_converged   = n(),
      n_extreme     = sum(extreme, na.rm = TRUE),
      prop_extreme  = n_extreme / n_converged,
      n_improper    = sum(improper, na.rm = TRUE),
      prop_improper = n_improper / n_converged,
      n_both        = sum(extreme & improper, na.rm = TRUE),
      n_excluded    = sum(excluded, na.rm = TRUE),
      bias_all      = mean(err, na.rm = TRUE),
      rmse_all      = sqrt(mean(err^2, na.rm = TRUE)),
      bias_excl     = mean(err[!excluded], na.rm = TRUE),
      rmse_excl     = sqrt(mean(err[!excluded]^2, na.rm = TRUE)),
      median_bias   = median(err, na.rm = TRUE),
      mdae          = median(abs(err), na.rm = TRUE),
      .groups = "drop")
  write.csv(cond_df, out("extreme_cond.csv"), row.names = FALSE)

  # ── 반복 단위 (대상 경계모수 중 하나라도 극단) ──
  rep_level <- ld %>%
    group_by(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    summarise(extreme_rep  = any(extreme, na.rm = TRUE),
              improper_rep = any(improper, na.rm = TRUE), .groups = "drop") %>%
    left_join(rep_df %>% select(cond_code, rep_id, starts_with("prop")), by = c("cond_code", "rep_id")) %>%
    left_join(ghost_info, by = c("cond_code", "rep_id"))

  cond_rep <- rep_level %>%
    group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
    summarise(n_converged = n(),
              n_extreme_rep = sum(extreme_rep), prop_extreme_rep = n_extreme_rep / n_converged,
              n_improper_rep = sum(improper_rep), prop_improper_rep = n_improper_rep / n_converged,
              .groups = "drop")
  write.csv(cond_rep, out("extreme_cond_rep.csv"), row.names = FALSE)

  # ── 원인 조합 ──
  cause_names <- c(c1_near_zero = "near_zero_dak", c2_sf_rev = "sf_reversal",
                   c3_sparse = "sparse_adjacent", c4_ghost = "ghost_fill", c5_se_fail = "se_fail")
  ext <- ld %>% filter(extreme)
  cause_tbl <- if (nrow(ext) == 0) data.frame() else {
    combo <- apply(ext[, names(cause_names)], 1, function(r) {
      x <- cause_names[as.logical(r)]; if (length(x) == 0) "none" else paste(x, collapse = "+")
    })
    as.data.frame(table(param = ext$param, causes = combo), stringsAsFactors = FALSE) %>%
      filter(Freq > 0) %>% group_by(param) %>%
      mutate(prop = Freq / sum(Freq)) %>% ungroup() %>% arrange(param, desc(Freq))
  }
  write.csv(cause_tbl, out("extreme_cause_table.csv"), row.names = FALSE)

  # ── 요약 ──
  sink(out("extreme_summary.txt"))
  tryCatch({
    cat(sprintf("NRM 극단 오차 진단 [%s] — 대상: %s\n\n", label, paste(params, collapse = ", ")))
    cat("판정 기준\n")
    cat(sprintf("  극단 오차   : |b_hat - b_true| >= %g\n", ERR_CUT))
    cat(sprintf("  부적절한 해 : dak <= 0 (채점함수 전치) 또는 |dak| < %g\n", DAK_CUT))
    cat(sprintf("  제외 지표   : %s 인 경계모수 추정값 제외\n",
                if (EXCLUDE_IMPROPER) "극단 오차 또는 부적절한 해" else "극단 오차"))
    cat("  근거: NRM 극단값 고정 기준 문헌 없음 → 부적절한 해 판정, 제외 수 보고,\n")
    cat("        중앙값 기반 지표 병기 (다분 IRT 시뮬레이션 관행)\n\n")

    cat(sprintf("수렴 반복 수: %d  |  경계모수 추정값 수: %d\n\n",
                n_distinct(paste(ld$cond_code, ld$rep_id)), nrow(ld)))
    cat("경계모수별 극단·부적절·제외 수\n")
    print(as.data.frame(ld %>% group_by(param) %>%
      summarise(n = n(), n_extreme = sum(extreme, na.rm = TRUE),
                prop_extreme = mean(extreme, na.rm = TRUE),
                n_improper = sum(improper, na.rm = TRUE),
                n_excluded = sum(excluded, na.rm = TRUE))), digits = 4)
    cat(sprintf("\n극단 반복 (대상 경계모수 중 하나라도): %d / %d (%.2f%%)\n\n",
                sum(rep_level$extreme_rep), nrow(rep_level), 100 * mean(rep_level$extreme_rep)))

    cat("극단 오차 × 부적절한 해 교차표 (경계모수 추정값 단위)\n")
    print(table(extreme = ld$extreme, improper = ld$improper, useNA = "ifany"))

    cat("\n원인 조합 (극단 오차 경계모수, 상위)\n")
    if (nrow(cause_tbl) > 0) print(as.data.frame(cause_tbl %>% group_by(param) %>% slice_head(n = 8)),
                                   row.names = FALSE, digits = 3)

    cat("\n유령 응답자 비교 (응답 파일을 읽은 반복만)\n")
    g <- rep_level %>% filter(!is.na(n_ghost))
    print(as.data.frame(g %>% group_by(extreme_rep) %>%
      summarise(n_reps = n(), prop_with_ghost = mean(n_ghost > 0), mean_n_ghost = mean(n_ghost))),
      digits = 4)

    cat("\n극단값이 Bias/RMSE에 준 영향 (경계모수별, 조건 평균)\n")
    print(as.data.frame(cond_df %>% group_by(param) %>%
      summarise(bias_all = mean(bias_all), bias_excl = mean(bias_excl, na.rm = TRUE),
                median_bias = mean(median_bias),
                rmse_all = mean(rmse_all), rmse_excl = mean(rmse_excl, na.rm = TRUE),
                mdae = mean(mdae))), digits = 4)
  }, finally = sink())

  # ── 그래프 1: 극단 반복 비율 (기존 전치 그래프 배치) ──
  p1 <- ggplot(add_factors(cond_rep), aes(x = iv3, y = prop_extreme_rep)) +
    geom_col(fill = "firebrick") +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "b-parameter Mean (IV3)",
         y = sprintf("Extreme-error rate (|error| >= %g, %s)", ERR_CUT, paste(params, collapse = ","))) +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out("extreme_rate_plot.png"), p1, width = 24, height = 18, dpi = 150)

  # ── 그래프 1-b: 경계모수별 극단 오차 비율 (전치 쌍별 막대 그래프와 같은 배치) ──
  param_colours <- c(b1 = "#4C9BE8", b2 = "#59A14F", b3 = "#F2A541", b4 = "#D1495B")
  p1b <- ggplot(add_factors(cond_df), aes(x = iv3, y = prop_extreme, fill = param)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_fill_manual(values = param_colours[params], name = "Boundary") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "b-parameter Mean (IV3)",
         y = sprintf("Extreme-error rate (|error| >= %g)", ERR_CUT)) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out("extreme_rate_by_param_bar.png"), p1b, width = 24, height = 18, dpi = 150)

  # ── 그래프 2: 범주 응답 비율 (극단 vs 정상 반복) ──
  p2_df <- rep_level %>% filter(!is.na(prop0)) %>%
    pivot_longer(prop0:prop4, names_to = "category", values_to = "prop") %>%
    mutate(category = factor(category, labels = paste("Cat", 1:5)),
           group = ifelse(extreme_rep, "Extreme", "Normal"))
  p2 <- ggplot(p2_df, aes(x = category, y = prop, fill = group)) +
    geom_boxplot(outlier.size = 0.4) +
    scale_fill_manual(values = c(Extreme = "firebrick", Normal = "grey70"), name = NULL) +
    labs(x = "Item 1 response category", y = "Proportion of responses (ghosts excluded)") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
  ggsave(out("extreme_resp_dist.png"), p2, width = 9, height = 6, dpi = 150)

  # ── 그래프 3: 양옆 범주 최소 응답 수 vs |오차| ──
  p3 <- ggplot(ld %>% filter(!is.na(n_min_adj), is.finite(err)),
               aes(x = n_min_adj + 1, y = abs(err) + 1e-3, colour = extreme)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_hline(yintercept = ERR_CUT, linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"), name = "Extreme") +
    facet_wrap(~ param) +
    labs(x = "min(adjacent category count) + 1  (log)", y = "|error|  (log)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(out("extreme_vs_sparse.png"), p3, width = 10, height = 7, dpi = 150)

  # ── 그래프 4: 대표 사례 범주 반응 곡선 ──
  cases <- ld %>% filter(extreme) %>% arrange(desc(abs(err))) %>%
    distinct(cond_code, rep_id, .keep_all = TRUE) %>% slice_head(n = N_CASES)
  theta <- seq(-4, 4, length.out = 401)
  for (i in seq_len(nrow(cases))) {
    cs <- cases[i, ]
    rr <- rep_df %>% filter(cond_code == cs$cond_code, rep_id == cs$rep_id)
    sf <- c(0, 1, 2, rr$iv1_sf4, 4)
    bt <- rr$iv3_b_mean + rr$iv2_b_interval * c(-1.5, -0.5, 0.5, 1.5)
    p_true <- gpcm_crc(theta, bt, sf)
    p_est  <- nrm_crc(theta, rr$a1, unlist(rr[paste0("ak", 0:4)]), unlist(rr[paste0("d", 0:4)]))
    cat_names <- paste("Cat", 1:5)
    p_true <- unname(p_true); colnames(p_true) <- cat_names
    p_est  <- unname(p_est);  colnames(p_est)  <- cat_names
    crc <- bind_rows(
      data.frame(theta = theta, p_true, model = "True (GPCM, generating)", check.names = FALSE),
      data.frame(theta = theta, p_est,  model = "Estimated (NRM)",         check.names = FALSE)) %>%
      pivot_longer(all_of(cat_names), names_to = "category", values_to = "p")
    sub <- sprintf("counts (ghosts excl.): %s | ghosts: %s (cats %s) | ak: %s | %s err = %.2f",
                   paste(unlist(rr[paste0("n", 0:4)]), collapse = "/"),
                   ifelse(is.na(cs$n_ghost), "NA", cs$n_ghost),
                   ifelse(is.na(cs$ghost_cats) | cs$ghost_cats == "", "-", cs$ghost_cats),
                   paste(sprintf("%.2f", unlist(rr[paste0("ak", 0:4)])), collapse = ","),
                   cs$param, cs$err)
    pc <- ggplot(crc, aes(theta, p, colour = category, linetype = model)) +
      geom_line(linewidth = 0.8) +
      geom_vline(xintercept = bt, colour = "grey40", linetype = "dotted") +
      labs(title = sprintf("cond %s rep %d", cs$cond_code, cs$rep_id), subtitle = sub,
           x = expression(theta), y = "P(category)", colour = NULL, linetype = NULL) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
    ggsave(out(sprintf("case_%s_rep%04d_%s.png", cs$cond_code, cs$rep_id, cs$param)),
           pc, width = 9, height = 6, dpi = 150)
  }

  # ── GLM: 극단 반복 여부 ~ IV1*IV2*IV3 ──
  sink(out("extreme_glm.txt"))
  tryCatch({
    gd <- add_factors(rep_level) %>% mutate(y = as.numeric(extreme_rep))
    ivs <- c("iv1", "iv2", "iv3")
    ivs <- ivs[sapply(ivs, function(v) nlevels(droplevels(gd[[v]])) >= 2)]
    cat(sprintf("로지스틱 GLM: 극단 반복 여부 (%s 중 |오차| >= %g) ~ %s\n",
                paste(params, collapse = ","), ERR_CUT, paste(ivs, collapse = " * ")))
    cat(sprintf("반복 수: %d | 극단 반복: %d (%.2f%%)\n\n", nrow(gd), sum(gd$y), 100 * mean(gd$y)))
    if (length(ivs) == 0 || var(gd$y) == 0) {
      cat("독립변수 수준 또는 결과 분산 부족 — GLM 생략\n")
    } else {
      fit <- suppressWarnings(glm(as.formula(paste("y ~", paste(ivs, collapse = " * "))),
                                  data = gd, family = binomial("logit")))
      lrt <- suppressWarnings(drop1(fit, scope = ~., test = "Chisq"))
      lrt_col <- intersect(c("LRT", "Chisq"), colnames(lrt))[1]
      tm  <- rownames(lrt)[rownames(lrt) != "<none>"]
      res <- as.data.frame(lrt[tm, , drop = FALSE])
      res$partial_pseudo_R2 <- res[[lrt_col]] / fit$null.deviance
      cat(sprintf("McFadden pseudo-R2: %.4f\n\n", 1 - deviance(fit) / fit$null.deviance))
      print(res, digits = 4)
    }
  }, error = function(e) cat("\n[오류]", conditionMessage(e), "\n"), finally = sink())

  message(sprintf("진단 [%s] 완료 → %s/", label, out_dir))
}

# =============================================================================
# 3. 실행
# =============================================================================
run_diag(c("b1", "b2", "b3", "b4"), "all_b")
run_diag(c("b3", "b4"),             "b34")

# =============================================================================
# 4. 전치 비율 분해 (output/nrm/diagnostics/all_b/)
#   전치 = 인접 경계모수 쌍 중 하나라도 b_k ≥ b_{k+1}  (GPCM·12번 NRM 전치 그래프와 같은 기준)
#   전치된 반복을 원인별로 나눔 — 분모: 수렴 반복 (b1~b4 모두 계산된 반복)
#     극단 오차    : b1~b4 중 하나라도 |오차| ≥ ERR_CUT
#     채점함수 이상: Δak1~Δak4 중 하나라도 Δak ≤ 0 (채점함수 전치) 또는 |Δak| < DAK_CUT (분모 붕괴)
#   범주: 둘 다 / 극단 오차만 / 채점함수 이상만 / 둘 다 아님(깨끗한 전치)
#   → 네 범주의 합 = 12번 nrm_transposition_plot.png 의 전치 비율
#   → 16번 excl_extreme 전치 비율       = 채점함수 이상만 + 둘 다 아님
#      16번 excl_extreme_sfrev 전치 비율 = 둘 다 아님
#
#   극단 오차는 대부분 분모 a1·Δak 가 0에 가까워서 생기므로 "극단 오차만" 범주는
#   실제 데이터에서 거의 없을 수 있음 → 범주별 반복 수를 transposition_decomposition_counts.txt 에 기록,
#   전치 반복 안에서의 원인 구성비(transposition_decomposition_share.png)도 따로 그림
# =============================================================================
dec_dir <- file.path(DIAG_DIR, "all_b")
dir.create(dec_dir, recursive = TRUE, showWarnings = FALSE)

DEC_LEVELS  <- c("Both (extreme + SF reversal/near-zero dak)", "Extreme error only",
                 "SF reversal/near-zero dak only", "Neither (clean reversal)")
dec_colours <- setNames(c("#6A3D9A", "#D1495B", "#F2A541", "#4C9BE8"), DEC_LEVELS)

dec_rep <- rep_df %>%
  filter(!is.na(b1), !is.na(b2), !is.na(b3), !is.na(b4)) %>%
  mutate(
    transposed = b1 >= b2 | b2 >= b3 | b3 >= b4,
    extreme    = pmax(abs(b1 - (iv3_b_mean + iv2_b_interval * -1.5)),
                      abs(b2 - (iv3_b_mean + iv2_b_interval * -0.5)),
                      abs(b3 - (iv3_b_mean + iv2_b_interval *  0.5)),
                      abs(b4 - (iv3_b_mean + iv2_b_interval *  1.5))) >= ERR_CUT,
    sf_rev     = {
      dak <- cbind(ak1 - ak0, ak2 - ak1, ak3 - ak2, ak4 - ak3)
      rowSums(dak <= 0 | abs(dak) < DAK_CUT) > 0
    },
    cause = case_when(
      !transposed         ~ NA_character_,
      extreme & sf_rev    ~ DEC_LEVELS[1],
      extreme & !sf_rev   ~ DEC_LEVELS[2],
      !extreme & sf_rev   ~ DEC_LEVELS[3],
      TRUE                ~ DEC_LEVELS[4])
  )

dec_cond <- dec_rep %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
  summarise(
    n_converged    = n(),
    n_transposed   = sum(transposed),
    n_both         = sum(cause %in% DEC_LEVELS[1]),
    n_extreme_only = sum(cause %in% DEC_LEVELS[2]),
    n_sfrev_only   = sum(cause %in% DEC_LEVELS[3]),
    n_clean        = sum(cause %in% DEC_LEVELS[4]),
    .groups = "drop") %>%
  mutate(prop_transposed = n_transposed / n_converged,
         prop_both = n_both / n_converged, prop_extreme_only = n_extreme_only / n_converged,
         prop_sfrev_only = n_sfrev_only / n_converged, prop_clean = n_clean / n_converged)
write.csv(dec_cond, file.path(dec_dir, "transposition_decomposition.csv"), row.names = FALSE)

# 범주별 반복 수 (그래프에 안 보이는 범주가 실제로 0인지 확인용)
sink(file.path(dec_dir, "transposition_decomposition_counts.txt"))
tryCatch({
  cat("NRM 전치(인접 쌍 중 하나라도) 원인 분해 — 전체 조건 합계\n")
  cat(sprintf("기준: 극단 오차 |오차| >= %g (b1~b4), 채점함수 이상 Δak <= 0 또는 |Δak| < %g (Δak1~Δak4)\n\n",
              ERR_CUT, DAK_CUT))
  tot <- colSums(dec_cond[, c("n_converged", "n_transposed", "n_both",
                              "n_extreme_only", "n_sfrev_only", "n_clean")])
  print(data.frame(n = tot, prop_of_converged = round(tot / tot[["n_converged"]], 5),
                   prop_of_transposed = round(tot / tot[["n_transposed"]], 5)))
  cat("\n전치되지 않은 반복의 극단 오차·채점함수 이상 (참고)\n")
  print(with(dec_rep[!dec_rep$transposed, ], table(extreme, sf_rev, useNA = "ifany")))
}, finally = sink())

dec_long <- dec_cond %>%
  select(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, n_transposed,
         n_both, n_extreme_only, n_sfrev_only, n_clean,
         prop_both, prop_extreme_only, prop_sfrev_only, prop_clean) %>%
  pivot_longer(c(n_both, n_extreme_only, n_sfrev_only, n_clean,
                 prop_both, prop_extreme_only, prop_sfrev_only, prop_clean),
               names_to = c(".value", "cause"), names_pattern = "^(n|prop)_(.*)$") %>%
  mutate(share = ifelse(n_transposed > 0, n / n_transposed, 0),
         cause = factor(cause, levels = c("both", "extreme_only", "sfrev_only", "clean"),
                        labels = DEC_LEVELS)) %>%
  add_factors()

dec_plot <- function(yvar, ylab, fname) {
  p <- ggplot(dec_long, aes(x = iv3, y = .data[[yvar]], fill = cause)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    facet_grid(rows = vars(iv1), cols = vars(iv2),
               labeller = labeller(iv1 = function(x) paste0("s4=", x),
                                   iv2 = function(x) paste0("Int=", x))) +
    scale_fill_manual(values = dec_colours, name = NULL, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "b-parameter Mean (IV3)", y = ylab) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 7), axis.text = element_text(size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(nrow = 1))
  ggsave(file.path(dec_dir, fname), p, width = 24, height = 18, dpi = 150)
}
dec_plot("prop",  "Reversal Rate (NRM, any adjacent pair), by cause",
         "transposition_decomposition.png")
dec_plot("share", "Share of reversed replications, by cause",
         "transposition_decomposition_share.png")
message(sprintf("저장 완료: %s/transposition_decomposition{,_share}.png / .csv / _counts.txt", dec_dir))
