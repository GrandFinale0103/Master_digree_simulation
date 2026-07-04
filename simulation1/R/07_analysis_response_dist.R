# =============================================================================
# 07_analysis_response_dist.R
# 응답 분포 및 조건별 전치 발생 모수 정리
#
# 출력 테이블 (3종):
#
#   [테이블 1] output/analysis/resp_dist_rep.csv
#     반복별 × 문항별 응답 범주 인원 수
#     (유령 응답자 제외, 각 문항에서 0~4 범주에 몇 명이 응답했는지)
#     열: cond_code, rep_id, IV 정보, item, n_persons,
#         n0~n4 (인원수), prop0~prop4 (비율)
#
#   [테이블 2] output/analysis/resp_dist_cond.csv
#     조건별 × 문항별 평균 응답 분포 (반복 간 평균)
#     열: cond_code, IV 정보, item, n_reps,
#         mean_n0~mean_n4 (평균 인원수), mean_prop0~mean_prop4 (평균 비율)
#
#   [테이블 3] output/analysis/transposition_pairs_cond.csv
#     조건별 전치 발생 경계모수 쌍 정리
#     (문항 1의 b1<b2, b2<b3, b3<b4 각 쌍의 전치 발생 횟수·비율)
#     열: cond_code, IV 정보, n_reps, n_converged, pct_converged,
#         n_trans_any, pct_trans_any,
#         n_trans_12, pct_trans_12,
#         n_trans_23, pct_trans_23,
#         n_trans_34, pct_trans_34
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

N_CAT   <- 5L   # 응답 범주 수 (0~4)
N_ITEMS <- 20L  # 문항 수

# ── 조건 코드 → 실제 값 매핑 (시뮬레이션 설계 변경 시 여기만 수정) ──────────
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")

decode_cond <- function(cond_code) {
  d <- strsplit(cond_code, "")[[1]]
  list(
    iv1_sf4        = IV1_MAP[d[1]],
    iv2_b_interval = IV2_MAP[d[2]],
    iv3_b_mean     = IV3_MAP[d[3]],
    iv4_theta_dist = IV4_MAP[d[4]]
  )
}

parse_fname <- function(path) {
  fname <- basename(path)
  cond  <- str_match(fname, "_cond([0-9]+)_")[, 2]
  rep   <- as.integer(str_match(fname, "_rep([0-9]+)_")[, 2])
  list(cond_code = cond, rep_id = rep)
}

# =============================================================================
# ■ 테이블 1 & 2: 응답 분포 (반복별, 조건별)
# =============================================================================

resp_files <- list.files(
  path       = "output/responses",
  pattern    = "_response\\.csv$",
  full.names = TRUE
)

if (length(resp_files) == 0) {
  warning("응답 파일이 없습니다 (output/responses/). 테이블 1·2는 생성하지 않습니다.")
} else {
  cat(sprintf("응답 파일 %d개 로딩 중 (테이블 1·2)...\n", length(resp_files)))

  read_resp_dist <- function(path) {
    meta <- parse_fname(path)
    iv   <- decode_cond(meta$cond_code)

    df <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) return(NULL)

    # 유령 응답자(theta = NA) 제외
    df <- df[!is.na(df$theta), , drop = FALSE]
    n_persons <- nrow(df)

    item_cols <- grep("^item", names(df), value = TRUE)

    rows <- lapply(seq_along(item_cols), function(i) {
      col    <- item_cols[i]
      vals   <- as.integer(df[[col]])
      counts <- tabulate(vals + 1L, nbins = N_CAT)   # 0→1, ..., 4→5번째 bin
      props  <- if (n_persons > 0) counts / n_persons else rep(NA_real_, N_CAT)

      data.frame(
        cond_code      = meta$cond_code,
        rep_id         = meta$rep_id,
        iv1_sf4        = iv$iv1_sf4,
        iv2_b_interval = iv$iv2_b_interval,
        iv3_b_mean     = iv$iv3_b_mean,
        iv4_theta_dist = iv$iv4_theta_dist,
        item           = i,
        n_persons      = n_persons,
        n0 = counts[1], n1 = counts[2], n2 = counts[3],
        n3 = counts[4], n4 = counts[5],
        prop0 = props[1], prop1 = props[2], prop2 = props[3],
        prop3 = props[4], prop4 = props[5],
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  rep_list <- lapply(resp_files, read_resp_dist)
  rep_df   <- do.call(rbind, rep_list[!sapply(rep_list, is.null)])
  rep_df   <- rep_df[order(rep_df$cond_code, rep_df$rep_id, rep_df$item), ]

  write.csv(rep_df, "output/analysis/resp_dist_rep.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/resp_dist_rep.csv (%d행)\n", nrow(rep_df)))

  # ── 테이블 2: 조건 × 문항별 평균 응답 분포 ──────────────────────────────────
  group_cols <- c("cond_code", "iv1_sf4", "iv2_b_interval",
                  "iv3_b_mean", "iv4_theta_dist", "item")

  cond_df <- rep_df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_reps        = n(),
      mean_n_persons = mean(n_persons, na.rm = TRUE),
      mean_n0 = mean(n0, na.rm = TRUE), mean_n1 = mean(n1, na.rm = TRUE),
      mean_n2 = mean(n2, na.rm = TRUE), mean_n3 = mean(n3, na.rm = TRUE),
      mean_n4 = mean(n4, na.rm = TRUE),
      mean_prop0 = mean(prop0, na.rm = TRUE), mean_prop1 = mean(prop1, na.rm = TRUE),
      mean_prop2 = mean(prop2, na.rm = TRUE), mean_prop3 = mean(prop3, na.rm = TRUE),
      mean_prop4 = mean(prop4, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cond_code, item)

  write.csv(cond_df, "output/analysis/resp_dist_cond.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/resp_dist_cond.csv (%d행)\n", nrow(cond_df)))
}

# =============================================================================
# ■ 테이블 3: 조건별 전치 발생 경계모수 쌍 정리
# =============================================================================

est_files <- list.files(
  path       = "output/estimated_params",
  pattern    = "_est_params\\.csv$",
  full.names = TRUE
)

if (length(est_files) == 0) {
  warning("추정 결과 파일이 없습니다 (output/estimated_params/). 테이블 3은 생성하지 않습니다.")
} else {
  cat(sprintf("\n추정 파라미터 파일 %d개 로딩 중 (테이블 3)...\n", length(est_files)))

  read_trans_pairs <- function(path) {
    meta <- parse_fname(path)
    iv   <- decode_cond(meta$cond_code)

    df <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (is.null(df)) return(NULL)

    # 추정 실패 파일 처리
    if ("success" %in% names(df) && !isTRUE(as.logical(df$success[1]))) {
      return(data.frame(
        cond_code      = meta$cond_code,
        rep_id         = meta$rep_id,
        iv1_sf4        = iv$iv1_sf4,
        iv2_b_interval = iv$iv2_b_interval,
        iv3_b_mean     = iv$iv3_b_mean,
        iv4_theta_dist = iv$iv4_theta_dist,
        converged      = FALSE,
        trans_any      = NA, trans_12 = NA, trans_23 = NA, trans_34 = NA,
        b1 = NA, b2 = NA, b3 = NA, b4 = NA,
        stringsAsFactors = FALSE
      ))
    }

    # 문항 1 행 추출
    item1 <- df[df$item %in% c("Item1", "item1", 1), ]
    if (nrow(item1) == 0) item1 <- df[1, ]

    b1 <- item1$b1[1]; b2 <- item1$b2[1]
    b3 <- item1$b3[1]; b4 <- item1$b4[1]

    converged <- as.logical(item1$converged[1])

    # 경계모수 쌍별 전치 여부 (b_k >= b_{k+1} 이면 전치)
    if (isTRUE(converged) && !anyNA(c(b1, b2, b3, b4))) {
      t12 <- as.integer(b1 >= b2)
      t23 <- as.integer(b2 >= b3)
      t34 <- as.integer(b3 >= b4)
      t_any <- as.integer(t12 | t23 | t34)
    } else {
      t12 <- t23 <- t34 <- t_any <- NA_integer_
    }

    data.frame(
      cond_code      = meta$cond_code,
      rep_id         = meta$rep_id,
      iv1_sf4        = iv$iv1_sf4,
      iv2_b_interval = iv$iv2_b_interval,
      iv3_b_mean     = iv$iv3_b_mean,
      iv4_theta_dist = iv$iv4_theta_dist,
      converged      = isTRUE(converged),
      trans_any      = t_any,
      trans_12       = t12,
      trans_23       = t23,
      trans_34       = t34,
      b1 = b1, b2 = b2, b3 = b3, b4 = b4,
      stringsAsFactors = FALSE
    )
  }

  trans_list <- lapply(est_files, read_trans_pairs)
  trans_raw  <- do.call(rbind, trans_list[!sapply(trans_list, is.null)])
  cat(sprintf("  총 %d행 로딩 완료\n", nrow(trans_raw)))

  # ── 조건별 집계 ──────────────────────────────────────────────────────────────
  group_cols2 <- c("cond_code", "iv1_sf4", "iv2_b_interval",
                   "iv3_b_mean", "iv4_theta_dist")

  pairs_cond <- trans_raw %>%
    group_by(across(all_of(group_cols2))) %>%
    summarise(
      n_reps        = n(),
      n_converged   = sum(converged, na.rm = TRUE),
      pct_converged = 100 * n_converged / n_reps,
      n_trans_any   = sum(trans_any,  na.rm = TRUE),
      pct_trans_any = 100 * n_trans_any / n_converged,
      n_trans_12    = sum(trans_12,   na.rm = TRUE),
      pct_trans_12  = 100 * n_trans_12 / n_converged,
      n_trans_23    = sum(trans_23,   na.rm = TRUE),
      pct_trans_23  = 100 * n_trans_23 / n_converged,
      n_trans_34    = sum(trans_34,   na.rm = TRUE),
      pct_trans_34  = 100 * n_trans_34 / n_converged,
      .groups = "drop"
    ) %>%
    arrange(cond_code)

  write.csv(pairs_cond, "output/analysis/transposition_pairs_cond.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/transposition_pairs_cond.csv (%d행)\n",
              nrow(pairs_cond)))

  # ── 반복별 원자료도 저장 (조건별 집계 전 원자료) ─────────────────────────────
  trans_raw_out <- trans_raw %>%
    select(cond_code, rep_id, iv1_sf4, iv2_b_interval, iv3_b_mean, iv4_theta_dist,
           converged, b1, b2, b3, b4,
           trans_any, trans_12, trans_23, trans_34) %>%
    arrange(cond_code, rep_id)

  write.csv(trans_raw_out, "output/analysis/transposition_pairs_rep.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/transposition_pairs_rep.csv (%d행)\n",
              nrow(trans_raw_out)))
}

# =============================================================================
# 최종 요약
# =============================================================================
cat("\n================================================================\n")
cat("전체 분석 완료\n")
cat("================================================================\n")

if (exists("rep_df")) {
  cat(sprintf("  [테이블 1] 반복별 응답 분포 : %d행 (반복 %d개 × 문항 %d개)\n",
              nrow(rep_df),
              n_distinct(rep_df[, c("cond_code", "rep_id")]),
              n_distinct(rep_df$item)))
  cat(sprintf("  [테이블 2] 조건별 평균 분포  : %d행 (조건 %d개 × 문항 %d개)\n",
              nrow(cond_df),
              n_distinct(cond_df$cond_code),
              n_distinct(cond_df$item)))
}
if (exists("pairs_cond")) {
  cat(sprintf("  [테이블 3] 조건별 전치 쌍    : %d행 (조건 %d개)\n",
              nrow(pairs_cond), nrow(pairs_cond)))
  cat(sprintf("             전체 전치 비율: %.1f%% ~ %.1f%%\n",
              min(pairs_cond$pct_trans_any, na.rm = TRUE),
              max(pairs_cond$pct_trans_any, na.rm = TRUE)))
  cat(sprintf("             b1≥b2: %.1f%% ~ %.1f%%  ",
              min(pairs_cond$pct_trans_12, na.rm = TRUE),
              max(pairs_cond$pct_trans_12, na.rm = TRUE)))
  cat(sprintf("b2≥b3: %.1f%% ~ %.1f%%  ",
              min(pairs_cond$pct_trans_23, na.rm = TRUE),
              max(pairs_cond$pct_trans_23, na.rm = TRUE)))
  cat(sprintf("b3≥b4: %.1f%% ~ %.1f%%\n",
              min(pairs_cond$pct_trans_34, na.rm = TRUE),
              max(pairs_cond$pct_trans_34, na.rm = TRUE)))
}

cat("\n출력 파일:\n")
cat("  output/analysis/resp_dist_rep.csv           (테이블 1)\n")
cat("  output/analysis/resp_dist_cond.csv          (테이블 2)\n")
cat("  output/analysis/transposition_pairs_rep.csv (테이블 3 원자료)\n")
cat("  output/analysis/transposition_pairs_cond.csv(테이블 3 집계)\n")
