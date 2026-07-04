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
# =============================================================================
est_files <- list.files(
  path = "output/estimated_params", pattern = "_est_params\\.csv$", full.names = TRUE
)

if (length(est_files) == 0) {
  warning("추정 결과 파일이 없습니다. 테이블 3은 생성하지 않습니다.")
} else {
  cat(sprintf("\n추정 파라미터 파일 %d개 로딩 중 (테이블 3)...\n", length(est_files)))

  USE_FREAD2 <- requireNamespace("data.table", quietly = TRUE)
  read_est   <- if (USE_FREAD2)
    function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  else
    function(p) read.csv(p, stringsAsFactors = FALSE)

  read_trans_pairs <- function(path) {
    meta <- parse_fname(path)
    iv   <- decode_cond(meta$cond_code)

    df <- tryCatch(read_est(path), error = function(e) NULL)
    if (is.null(df)) return(NULL)

    if ("success" %in% names(df) && !isTRUE(as.logical(df$success[1]))) {
      return(data.frame(
        cond_code=meta$cond_code, rep_id=meta$rep_id,
        iv1_sf4=iv$iv1_sf4, iv2_b_interval=iv$iv2_b_interval,
        iv3_b_mean=iv$iv3_b_mean, iv4_theta_dist=iv$iv4_theta_dist,
        converged=FALSE,
        b1=NA_real_, b2=NA_real_, b3=NA_real_, b4=NA_real_,
        trans_any=NA_integer_, trans_12=NA_integer_,
        trans_23=NA_integer_,  trans_34=NA_integer_,
        stringsAsFactors=FALSE
      ))
    }

    item1 <- df[df$item %in% c("Item1","item1",1), ]
    if (nrow(item1) == 0) item1 <- df[1, ]

    b1 <- item1$b1[1]; b2 <- item1$b2[1]
    b3 <- item1$b3[1]; b4 <- item1$b4[1]
    conv <- as.logical(item1$converged[1])

    if (isTRUE(conv) && !anyNA(c(b1,b2,b3,b4))) {
      t12 <- as.integer(b1 >= b2); t23 <- as.integer(b2 >= b3)
      t34 <- as.integer(b3 >= b4); tany <- as.integer(t12|t23|t34)
    } else {
      t12 <- t23 <- t34 <- tany <- NA_integer_
    }

    data.frame(
      cond_code=meta$cond_code, rep_id=meta$rep_id,
      iv1_sf4=iv$iv1_sf4, iv2_b_interval=iv$iv2_b_interval,
      iv3_b_mean=iv$iv3_b_mean, iv4_theta_dist=iv$iv4_theta_dist,
      converged=isTRUE(conv),
      b1=b1, b2=b2, b3=b3, b4=b4,
      trans_any=tany, trans_12=t12, trans_23=t23, trans_34=t34,
      stringsAsFactors=FALSE
    )
  }

  # est_params 파일은 작으므로 단순 lapply 사용
  trans_list <- lapply(est_files, read_trans_pairs)
  trans_raw  <- do.call(rbind, trans_list[!sapply(trans_list, is.null)])
  cat(sprintf("  총 %d행 로딩 완료\n", nrow(trans_raw)))

  trans_raw <- trans_raw[order(trans_raw$cond_code, trans_raw$rep_id), ]
  write.csv(trans_raw, "output/analysis/transposition_pairs_rep.csv", row.names = FALSE)
  cat(sprintf("  저장 완료: output/analysis/transposition_pairs_rep.csv (%d행)\n",
              nrow(trans_raw)))

  grp <- c("cond_code","iv1_sf4","iv2_b_interval","iv3_b_mean","iv4_theta_dist")
  pairs_cond <- do.call(rbind, lapply(split(trans_raw, trans_raw$cond_code), function(d) {
    nc  <- sum(d$converged, na.rm=TRUE)
    data.frame(
      cond_code      = d$cond_code[1],
      iv1_sf4        = d$iv1_sf4[1],
      iv2_b_interval = d$iv2_b_interval[1],
      iv3_b_mean     = d$iv3_b_mean[1],
      iv4_theta_dist = d$iv4_theta_dist[1],
      n_reps         = nrow(d),
      n_converged    = nc,
      pct_converged  = 100*nc/nrow(d),
      n_trans_any    = sum(d$trans_any,  na.rm=TRUE),
      pct_trans_any  = 100*sum(d$trans_any,  na.rm=TRUE)/nc,
      n_trans_12     = sum(d$trans_12,   na.rm=TRUE),
      pct_trans_12   = 100*sum(d$trans_12,   na.rm=TRUE)/nc,
      n_trans_23     = sum(d$trans_23,   na.rm=TRUE),
      pct_trans_23   = 100*sum(d$trans_23,   na.rm=TRUE)/nc,
      n_trans_34     = sum(d$trans_34,   na.rm=TRUE),
      pct_trans_34   = 100*sum(d$trans_34,   na.rm=TRUE)/nc,
      stringsAsFactors = FALSE
    )
  }))
  pairs_cond <- pairs_cond[order(pairs_cond$cond_code), ]
  write.csv(pairs_cond, "output/analysis/transposition_pairs_cond.csv", row.names=FALSE)
  cat(sprintf("  저장 완료: output/analysis/transposition_pairs_cond.csv (%d행)\n",
              nrow(pairs_cond)))
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
