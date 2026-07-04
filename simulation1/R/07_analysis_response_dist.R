# =============================================================================
# 07_analysis_response_dist.R
# 응답 분포 및 조건별 전치 발생 모수 정리
#
# 성능 최적화:
#   - data.table::fread 사용 (read.csv 대비 5~10배 빠름, 없으면 자동 폴백)
#   - parallel 패키지로 멀티코어 병렬 처리
#   - Table 2: 누적 합계(accumulator) 패턴 — 원자료를 메모리에 쌓지 않음
#   - Table 1: 배치(BATCH_SIZE 파일) 단위 append 쓰기 — GC 주기적 실행
#
# 출력 테이블 (4종):
#   [테이블 1] output/analysis/resp_dist_rep.csv
#     반복별 × 문항별 응답 범주 인원 수 (유령 응답자 제외)
#     열: cond_code, rep_id, IV 정보, item, n_persons, n0~n4, prop0~prop4
#
#   [테이블 2] output/analysis/resp_dist_cond.csv
#     조건별 × 문항별 평균 응답 분포 (반복 간 평균)
#     열: cond_code, IV 정보, item, n_reps, mean_n_persons,
#         mean_n0~mean_n4, mean_prop0~mean_prop4
#
#   [테이블 3a] output/analysis/transposition_pairs_rep.csv
#     반복별 경계모수 쌍 전치 여부 원자료
#   [테이블 3b] output/analysis/transposition_pairs_cond.csv
#     조건별 경계모수 쌍 전치 발생 집계
# =============================================================================

suppressPackageStartupMessages(library(parallel))

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── 성능 파라미터 (필요 시 조정) ─────────────────────────────────────────────
BATCH_SIZE <- 200L   # 배치당 처리 파일 수 (메모리 ↔ 속도 균형)
N_CORES    <- max(1L, detectCores(logical = FALSE) - 1L)  # 물리 코어 기준

# ── 설계 상수 ─────────────────────────────────────────────────────────────────
N_CAT   <- 5L
N_ITEMS <- 20L

IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")

decode_cond <- function(cond_code) {
  d <- strsplit(cond_code, "")[[1]]
  list(
    iv1_sf4        = unname(IV1_MAP[d[1]]),
    iv2_b_interval = unname(IV2_MAP[d[2]]),
    iv3_b_mean     = unname(IV3_MAP[d[3]]),
    iv4_theta_dist = unname(IV4_MAP[d[4]])
  )
}

parse_fname <- function(path) {
  b <- basename(path)
  list(
    cond_code = sub(".*_cond([0-9]+)_.*", "\\1", b),
    rep_id    = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", b))
  )
}

# =============================================================================
# ■ 테이블 1 & 2: 응답 분포
# =============================================================================
resp_files <- list.files(
  path = "output/responses", pattern = "_response\\.csv$", full.names = TRUE
)

if (length(resp_files) == 0) {
  warning("응답 파일이 없습니다 (output/responses/). 테이블 1·2는 생성하지 않습니다.")
} else {

  # ── CSV 고속 읽기: data.table::fread 우선, 없으면 read.csv ──────────────────
  USE_FREAD <- requireNamespace("data.table", quietly = TRUE)
  if (USE_FREAD) {
    cat("  data.table::fread 사용 (고속 모드)\n")
    read_fast <- function(p)
      data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  } else {
    cat("  read.csv 사용 (data.table 없음 — install.packages('data.table') 권장)\n")
    read_fast <- function(p) read.csv(p, stringsAsFactors = FALSE)
  }

  n_files <- length(resp_files)
  cat(sprintf("응답 파일 %d개 처리 시작 (코어 %d개, 배치 %d파일)\n",
              n_files, N_CORES, BATCH_SIZE))

  # ── Table 1: 헤더만 먼저 파일에 쓰기 (이후 append) ─────────────────────────
  T1_PATH <- "output/analysis/resp_dist_rep.csv"
  T1_COLS <- c("cond_code", "rep_id", "iv1_sf4", "iv2_b_interval",
               "iv3_b_mean", "iv4_theta_dist", "item", "n_persons",
               paste0("n",    0:4),
               paste0("prop", 0:4))
  writeLines(paste(T1_COLS, collapse = ","), T1_PATH)

  # ── Table 2: 누적 합계 accumulator ──────────────────────────────────────────
  # key = "cond_code|item" → list(meta, sum_n[5], sum_prop[5], sum_persons, count)
  acc_env <- new.env(hash = TRUE, parent = emptyenv())

  # ── 단일 파일 처리 함수 (worker 내부·직접 호출 공용) ─────────────────────────
  process_one_resp <- function(path, read_fn, n_cat, iv1, iv2, iv3, iv4,
                               decode_fn, parse_fn) {
    meta <- parse_fn(path)
    iv   <- decode_fn(meta$cond_code)

    df <- tryCatch(read_fn(path), error = function(e) NULL)
    if (is.null(df)) return(NULL)

    df <- df[!is.na(df$theta), , drop = FALSE]
    n_persons <- nrow(df)
    if (n_persons == 0L) return(NULL)

    item_cols <- grep("^item", names(df), value = TRUE)
    n_it      <- length(item_cols)

    # 전 문항 벡터화 처리: [n_cat × n_it] 행렬
    counts_mat <- vapply(df[item_cols], function(col)
      tabulate(as.integer(col) + 1L, nbins = n_cat), integer(n_cat))
    props_mat  <- counts_mat / n_persons

    # Table 1 데이터: character 행렬 (write.table 속도 최적화)
    t1_mat <- cbind(
      cond_code      = meta$cond_code,
      rep_id         = meta$rep_id,
      iv1_sf4        = iv$iv1_sf4,
      iv2_b_interval = iv$iv2_b_interval,
      iv3_b_mean     = iv$iv3_b_mean,
      iv4_theta_dist = iv$iv4_theta_dist,
      item           = seq_len(n_it),
      n_persons      = n_persons,
      t(counts_mat),           # n0..n4
      t(round(props_mat, 6))   # prop0..prop4
    )

    # Table 2 accumulator 갱신 데이터 (리스트로 반환)
    acc_updates <- lapply(seq_len(n_it), function(j)
      list(
        key            = paste0(meta$cond_code, "|", j),
        cond_code      = meta$cond_code,
        item           = j,
        iv1_sf4        = iv$iv1_sf4,
        iv2_b_interval = iv$iv2_b_interval,
        iv3_b_mean     = iv$iv3_b_mean,
        iv4_theta_dist = iv$iv4_theta_dist,
        counts         = counts_mat[, j],
        props          = props_mat[, j],
        n_persons      = n_persons
      )
    )

    list(t1_mat = t1_mat, acc_updates = acc_updates)
  }

  # ── 배치 처리 루프 ────────────────────────────────────────────────────────────
  n_batches    <- ceiling(n_files / BATCH_SIZE)
  total_t1_rows <- 0L

  for (b in seq_len(n_batches)) {
    idx_s  <- (b - 1L) * BATCH_SIZE + 1L
    idx_e  <- min(b * BATCH_SIZE, n_files)
    batch  <- resp_files[idx_s:idx_e]

    # 병렬 실행 (코어 1개면 lapply로 폴백)
    if (N_CORES > 1L) {
      cl <- makeCluster(N_CORES)
      on.exit(stopCluster(cl), add = TRUE)
      clusterExport(cl,
        c("read_fast", "USE_FREAD", "N_CAT", "IV1_MAP", "IV2_MAP",
          "IV3_MAP", "IV4_MAP", "decode_cond", "parse_fname",
          "process_one_resp"),
        envir = environment())
      results <- parLapply(cl, batch, function(p)
        process_one_resp(p, read_fast, N_CAT,
                         IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                         decode_cond, parse_fname))
      stopCluster(cl)
      on.exit(NULL)  # 이미 실행됐으므로 해제
    } else {
      results <- lapply(batch, function(p)
        process_one_resp(p, read_fast, N_CAT,
                         IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                         decode_cond, parse_fname))
    }

    results <- results[!sapply(results, is.null)]
    if (length(results) == 0L) next

    # ── Table 1: 배치 append 쓰기 ───────────────────────────────────────────
    t1_batch <- do.call(rbind, lapply(results, `[[`, "t1_mat"))
    write.table(t1_batch, T1_PATH,
                append    = TRUE,
                sep       = ",",
                row.names = FALSE,
                col.names = FALSE,
                quote     = FALSE)
    total_t1_rows <- total_t1_rows + nrow(t1_batch)

    # ── Table 2 accumulator 갱신 ────────────────────────────────────────────
    for (res in results) {
      for (upd in res$acc_updates) {
        key <- upd$key
        if (exists(key, envir = acc_env, inherits = FALSE)) {
          old <- get(key, envir = acc_env)
          old$sum_n       <- old$sum_n       + upd$counts
          old$sum_prop    <- old$sum_prop    + upd$props
          old$sum_persons <- old$sum_persons + upd$n_persons
          old$count       <- old$count       + 1L
          assign(key, old, envir = acc_env)
        } else {
          assign(key, list(
            cond_code      = upd$cond_code,
            item           = upd$item,
            iv1_sf4        = upd$iv1_sf4,
            iv2_b_interval = upd$iv2_b_interval,
            iv3_b_mean     = upd$iv3_b_mean,
            iv4_theta_dist = upd$iv4_theta_dist,
            sum_n          = upd$counts,
            sum_prop       = upd$props,
            sum_persons    = upd$n_persons,
            count          = 1L
          ), envir = acc_env)
        }
      }
    }

    rm(results, t1_batch)
    gc(verbose = FALSE)

    cat(sprintf("\r  진행: %d / %d 배치 완료 (%d / %d 파일)...",
                b, n_batches, idx_e, n_files))
  }
  cat("\n")

  cat(sprintf("  저장 완료: %s (%d행)\n", T1_PATH, total_t1_rows))

  # ── Table 2: accumulator → 평균 계산 후 저장 ─────────────────────────────
  acc_keys <- ls(envir = acc_env)
  cond_df  <- do.call(rbind, lapply(acc_keys, function(key) {
    a <- get(key, envir = acc_env)
    data.frame(
      cond_code      = a$cond_code,
      item           = a$item,
      iv1_sf4        = a$iv1_sf4,
      iv2_b_interval = a$iv2_b_interval,
      iv3_b_mean     = a$iv3_b_mean,
      iv4_theta_dist = a$iv4_theta_dist,
      n_reps         = a$count,
      mean_n_persons = a$sum_persons / a$count,
      mean_n0  = a$sum_n[1] / a$count,
      mean_n1  = a$sum_n[2] / a$count,
      mean_n2  = a$sum_n[3] / a$count,
      mean_n3  = a$sum_n[4] / a$count,
      mean_n4  = a$sum_n[5] / a$count,
      mean_prop0 = a$sum_prop[1] / a$count,
      mean_prop1 = a$sum_prop[2] / a$count,
      mean_prop2 = a$sum_prop[3] / a$count,
      mean_prop3 = a$sum_prop[4] / a$count,
      mean_prop4 = a$sum_prop[5] / a$count,
      stringsAsFactors = FALSE
    )
  }))

  cond_df <- cond_df[order(cond_df$cond_code, cond_df$item), ]
  write.csv(cond_df, "output/analysis/resp_dist_cond.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/resp_dist_cond.csv (%d행)\n", nrow(cond_df)))

  rm(acc_env); gc(verbose = FALSE)
}

# =============================================================================
# ■ 테이블 3: 조건별 전치 발생 경계모수 쌍 정리 (est_params 파일 기반)
#   응답 파일과 동일한 배치-병렬-accumulator 패턴 적용
# =============================================================================
est_files <- list.files(
  path = "output/estimated_params", pattern = "_est_params\\.csv$", full.names = TRUE
)

if (length(est_files) == 0) {
  warning("추정 결과 파일이 없습니다. 테이블 3은 생성하지 않습니다.")
} else {

  USE_FREAD3 <- requireNamespace("data.table", quietly = TRUE)
  read_est <- if (USE_FREAD3)
    function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  else
    function(p) read.csv(p, stringsAsFactors = FALSE)

  n_est     <- length(est_files)
  n_bat_est <- ceiling(n_est / BATCH_SIZE)
  cat(sprintf("\n추정 파라미터 파일 %d개 처리 시작 (코어 %d개, 배치 %d파일)\n",
              n_est, N_CORES, BATCH_SIZE))

  # ── 단일 est_params 파일 처리 함수 ───────────────────────────────────────────
  process_one_est <- function(path, read_fn, iv1, iv2, iv3, iv4,
                              decode_fn, parse_fn) {
    meta <- parse_fn(path)
    iv   <- decode_fn(meta$cond_code)

    df <- tryCatch(read_fn(path), error = function(e) NULL)
    if (is.null(df)) return(NULL)

    # 추정 실패 파일 처리
    if ("success" %in% names(df) && !isTRUE(as.logical(df$success[1]))) {
      return(list(
        row = c(meta$cond_code, meta$rep_id,
                iv$iv1_sf4, iv$iv2_b_interval, iv$iv3_b_mean, iv$iv4_theta_dist,
                "FALSE", NA,NA,NA,NA, NA,NA,NA,NA),
        cond_code = meta$cond_code, iv = iv,
        conv = FALSE, tany = NA, t12 = NA, t23 = NA, t34 = NA
      ))
    }

    item1 <- df[df$item %in% c("Item1","item1",1), ]
    if (nrow(item1) == 0L) item1 <- df[1L, ]

    b1 <- item1$b1[1]; b2 <- item1$b2[1]
    b3 <- item1$b3[1]; b4 <- item1$b4[1]
    conv <- isTRUE(as.logical(item1$converged[1]))

    if (conv && !anyNA(c(b1, b2, b3, b4))) {
      t12 <- as.integer(b1 >= b2); t23 <- as.integer(b2 >= b3)
      t34 <- as.integer(b3 >= b4); tany <- as.integer(t12 | t23 | t34)
    } else {
      t12 <- t23 <- t34 <- tany <- NA_integer_
    }

    list(
      # Table 3a 행: character vector (write.table 최적화)
      row = c(meta$cond_code, meta$rep_id,
              iv$iv1_sf4, iv$iv2_b_interval, iv$iv3_b_mean, iv$iv4_theta_dist,
              as.character(conv), b1, b2, b3, b4, tany, t12, t23, t34),
      # Table 3b accumulator 갱신 정보
      cond_code = meta$cond_code, iv = iv,
      conv = conv, tany = tany, t12 = t12, t23 = t23, t34 = t34
    )
  }

  # ── Table 3a: 헤더 먼저 쓰기 ─────────────────────────────────────────────────
  T3_PATH <- "output/analysis/transposition_pairs_rep.csv"
  T3_COLS <- c("cond_code","rep_id","iv1_sf4","iv2_b_interval",
               "iv3_b_mean","iv4_theta_dist","converged",
               "b1","b2","b3","b4","trans_any","trans_12","trans_23","trans_34")
  writeLines(paste(T3_COLS, collapse = ","), T3_PATH)

  # ── Table 3b: 조건별 누적 합계 accumulator ───────────────────────────────────
  # key = cond_code → list(iv 정보, n_reps, n_conv, n_any, n_12, n_23, n_34)
  acc3_env <- new.env(hash = TRUE, parent = emptyenv())

  total_t3_rows <- 0L

  for (b in seq_len(n_bat_est)) {
    idx_s <- (b - 1L) * BATCH_SIZE + 1L
    idx_e <- min(b * BATCH_SIZE, n_est)
    batch <- est_files[idx_s:idx_e]

    # 병렬 실행
    if (N_CORES > 1L) {
      cl <- makeCluster(N_CORES)
      on.exit(stopCluster(cl), add = TRUE)
      clusterExport(cl,
        c("read_est", "USE_FREAD3", "IV1_MAP", "IV2_MAP", "IV3_MAP", "IV4_MAP",
          "decode_cond", "parse_fname", "process_one_est"),
        envir = environment())
      results <- parLapply(cl, batch, function(p)
        process_one_est(p, read_est, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                        decode_cond, parse_fname))
      stopCluster(cl)
      on.exit(NULL)
    } else {
      results <- lapply(batch, function(p)
        process_one_est(p, read_est, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                        decode_cond, parse_fname))
    }

    results <- results[!sapply(results, is.null)]
    if (length(results) == 0L) next

    # ── Table 3a: 배치 append 쓰기 ─────────────────────────────────────────
    t3_batch <- do.call(rbind, lapply(results, `[[`, "row"))
    write.table(t3_batch, T3_PATH,
                append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    total_t3_rows <- total_t3_rows + nrow(t3_batch)

    # ── Table 3b accumulator 갱신 ───────────────────────────────────────────
    for (res in results) {
      key <- res$cond_code
      if (exists(key, envir = acc3_env, inherits = FALSE)) {
        old <- get(key, envir = acc3_env)
        old$n_reps <- old$n_reps + 1L
        old$n_conv <- old$n_conv + as.integer(isTRUE(res$conv))
        old$n_any  <- old$n_any  + ifelse(is.na(res$tany), 0L, res$tany)
        old$n_12   <- old$n_12   + ifelse(is.na(res$t12),  0L, res$t12)
        old$n_23   <- old$n_23   + ifelse(is.na(res$t23),  0L, res$t23)
        old$n_34   <- old$n_34   + ifelse(is.na(res$t34),  0L, res$t34)
        assign(key, old, envir = acc3_env)
      } else {
        assign(key, list(
          cond_code      = res$cond_code,
          iv1_sf4        = res$iv$iv1_sf4,
          iv2_b_interval = res$iv$iv2_b_interval,
          iv3_b_mean     = res$iv$iv3_b_mean,
          iv4_theta_dist = res$iv$iv4_theta_dist,
          n_reps = 1L,
          n_conv = as.integer(isTRUE(res$conv)),
          n_any  = ifelse(is.na(res$tany), 0L, res$tany),
          n_12   = ifelse(is.na(res$t12),  0L, res$t12),
          n_23   = ifelse(is.na(res$t23),  0L, res$t23),
          n_34   = ifelse(is.na(res$t34),  0L, res$t34)
        ), envir = acc3_env)
      }
    }

    rm(results, t3_batch)
    gc(verbose = FALSE)

    cat(sprintf("\r  진행: %d / %d 배치 완료 (%d / %d 파일)...",
                b, n_bat_est, idx_e, n_est))
  }
  cat("\n")
  cat(sprintf("  저장 완료: %s (%d행)\n", T3_PATH, total_t3_rows))

  # ── Table 3b: accumulator → 집계 후 저장 ─────────────────────────────────
  pairs_cond <- do.call(rbind, lapply(ls(envir = acc3_env), function(key) {
    a  <- get(key, envir = acc3_env)
    nc <- a$n_conv
    data.frame(
      cond_code      = a$cond_code,
      iv1_sf4        = a$iv1_sf4,
      iv2_b_interval = a$iv2_b_interval,
      iv3_b_mean     = a$iv3_b_mean,
      iv4_theta_dist = a$iv4_theta_dist,
      n_reps         = a$n_reps,
      n_converged    = nc,
      pct_converged  = 100 * nc / a$n_reps,
      n_trans_any    = a$n_any,
      pct_trans_any  = if (nc > 0) 100 * a$n_any / nc else NA_real_,
      n_trans_12     = a$n_12,
      pct_trans_12   = if (nc > 0) 100 * a$n_12  / nc else NA_real_,
      n_trans_23     = a$n_23,
      pct_trans_23   = if (nc > 0) 100 * a$n_23  / nc else NA_real_,
      n_trans_34     = a$n_34,
      pct_trans_34   = if (nc > 0) 100 * a$n_34  / nc else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  pairs_cond <- pairs_cond[order(pairs_cond$cond_code), ]
  write.csv(pairs_cond, "output/analysis/transposition_pairs_cond.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/transposition_pairs_cond.csv (%d행)\n",
              nrow(pairs_cond)))

  rm(acc3_env); gc(verbose = FALSE)
}

# =============================================================================
# 최종 요약
# =============================================================================
cat("\n================================================================\n")
cat("전체 분석 완료\n")
cat("================================================================\n")
if (exists("total_t1_rows")) {
  cat(sprintf("  [테이블 1] 반복별 응답 분포 : %d행\n", total_t1_rows))
  cat(sprintf("  [테이블 2] 조건별 평균 분포 : %d행\n", nrow(cond_df)))
}
if (exists("pairs_cond")) {
  cat(sprintf("  [테이블 3] 조건별 전치 쌍   : %d행\n", nrow(pairs_cond)))
  cat(sprintf("             전치 비율 범위: %.1f%% ~ %.1f%%\n",
              min(pairs_cond$pct_trans_any, na.rm=TRUE),
              max(pairs_cond$pct_trans_any, na.rm=TRUE)))
}
cat("\n출력 파일:\n")
cat("  output/analysis/resp_dist_rep.csv\n")
cat("  output/analysis/resp_dist_cond.csv\n")
cat("  output/analysis/transposition_pairs_rep.csv\n")
cat("  output/analysis/transposition_pairs_cond.csv\n")
