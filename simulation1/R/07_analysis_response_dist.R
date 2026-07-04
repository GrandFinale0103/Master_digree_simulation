# =============================================================================
# 07_analysis_response_dist.R
# 응답 분포 및 조건별 전치 발생 모수 정리
#
# [파일 검증]
#   - 기대 조건 코드(108개) 대비 실제 파일 존재 여부를 검사
#   - 누락 조건, 반복 횟수 미달 조건, 응시자 수 분포를 로그에 기록
#   - 존재하는 파일만으로 유연하게 집계 (N_PERSONS·N_REPS 고정 가정 없음)
#
# [집계 범위]
#   - 테이블 1·2: 문항 1(조작 문항)만 집계
#   - 테이블 3: 문항 1 경계모수 쌍(b1<b2, b2<b3, b3<b4) 전치 여부
#
# [성능]
#   - data.table::fread 우선 사용 (없으면 read.csv)
#   - parallel::parLapply 병렬 처리
#   - 배치 append 쓰기 + accumulator 패턴으로 메모리 절약
#
# 출력 파일:
#   output/analysis/07_log.txt                    — 실행 로그
#   output/analysis/resp_dist_rep_item1.csv       — 문항 1 반복별 응답 분포
#   output/analysis/resp_dist_cond_item1.csv      — 문항 1 조건별 평균 응답 분포
#   output/analysis/transposition_pairs_rep.csv   — 반복별 전치 쌍 원자료
#   output/analysis/transposition_pairs_cond.csv  — 조건별 전치 쌍 집계
# =============================================================================

suppressPackageStartupMessages(library(parallel))

dir.create("output/analysis", recursive = TRUE, showWarnings = FALSE)

# ── 성능 파라미터 ─────────────────────────────────────────────────────────────
BATCH_SIZE <- 200L
N_CORES    <- max(1L, detectCores(logical = FALSE) - 1L)

# ── 설계 상수 (IV 매핑) ───────────────────────────────────────────────────────
N_CAT   <- 5L    # 응답 범주 수 (0~4)
IV1_MAP <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
IV4_MAP <- c("1" = "pos_skew", "2" = "normal", "3" = "neg_skew", "4" = "uniform")

# 기대 조건 코드 전체 목록 (3×3×3×4 = 108개)
EXPECTED_CONDS <- character(108)
k <- 0L
for (d1 in 1:3) for (d2 in 1:3) for (d3 in 1:3) for (d4 in 1:4) {
  k <- k + 1L; EXPECTED_CONDS[k] <- paste0(d1, d2, d3, d4)
}

# ── 로그 시스템 ───────────────────────────────────────────────────────────────
LOG_PATH  <- "output/analysis/07_log.txt"
log_lines <- character(0)

lg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}
lg_section <- function(title) lg(sprintf("\n[%s]  %s", Sys.time(), title))
flush_log  <- function() writeLines(log_lines, LOG_PATH)

lg("================================================================")
lg(sprintf("07_analysis_response_dist.R  실행 시작: %s", Sys.time()))
lg(sprintf("작업 디렉토리: %s", getwd()))
lg("================================================================")

# ── 헬퍼 함수 ─────────────────────────────────────────────────────────────────
parse_fname <- function(path) {
  b <- basename(path)
  list(
    cond_code = sub(".*_cond([0-9]+)_.*", "\\1", b),
    rep_id    = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", b))
  )
}

decode_cond <- function(cond_code) {
  d <- strsplit(cond_code, "")[[1]]
  list(
    iv1_sf4        = unname(IV1_MAP[d[1]]),
    iv2_b_interval = unname(IV2_MAP[d[2]]),
    iv3_b_mean     = unname(IV3_MAP[d[3]]),
    iv4_theta_dist = unname(IV4_MAP[d[4]])
  )
}

# ── 공통 파일 검증 함수 ───────────────────────────────────────────────────────
validate_files <- function(files, section_label) {
  lg_section(paste(section_label, "파일 검증"))

  if (length(files) == 0) {
    lg("  파일 없음 — 이 섹션을 건너뜁니다.")
    return(invisible(NULL))
  }

  found_conds <- sapply(files, function(p) parse_fname(p)$cond_code)
  reps_tbl    <- sort(table(found_conds))          # 조건별 파일(반복) 수
  uniq_conds  <- names(reps_tbl)

  missing_conds  <- setdiff(EXPECTED_CONDS, uniq_conds)
  extra_conds    <- setdiff(uniq_conds, EXPECTED_CONDS)   # 예상 외 조건코드
  max_reps       <- max(reps_tbl)
  under_conds    <- reps_tbl[reps_tbl < max_reps]         # 반복 횟수 미달

  lg(sprintf("  전체 파일 수      : %d개", length(files)))
  lg(sprintf("  발견된 조건 수    : %d / 108개", length(uniq_conds)))
  lg(sprintf("  누락 조건 수      : %d개", length(missing_conds)))
  lg(sprintf("  예상 외 조건 수   : %d개", length(extra_conds)))
  lg(sprintf("  최다 반복 횟수    : %d회 (기준)", max_reps))
  lg(sprintf("  반복 미달 조건 수 : %d개", length(under_conds)))

  if (length(missing_conds) > 0) {
    lg("  ※ 누락 조건 코드:")
    for (chunk in split(missing_conds, ceiling(seq_along(missing_conds) / 18)))
      lg("    ", paste(chunk, collapse = " "))
  }
  if (length(extra_conds) > 0) {
    lg("  ※ 예상 외 조건 코드: ", paste(extra_conds, collapse = " "))
  }
  if (length(under_conds) > 0) {
    lg("  ※ 반복 미달 조건 (조건코드=횟수):")
    under_str <- paste(sprintf("%s=%d", names(under_conds), as.integer(under_conds)),
                       collapse = "  ")
    lg("    ", under_str)
  }

  flush_log()
  invisible(list(
    reps_tbl      = reps_tbl,
    missing_conds = missing_conds,
    extra_conds   = extra_conds,
    under_conds   = under_conds,
    max_reps      = max_reps
  ))
}

# ── fread / read.csv 선택 ─────────────────────────────────────────────────────
USE_FREAD <- requireNamespace("data.table", quietly = TRUE)
read_fast <- if (USE_FREAD) {
  function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
} else {
  function(p) read.csv(p, stringsAsFactors = FALSE)
}
lg(sprintf("\nCSV 읽기: %s", if (USE_FREAD) "data.table::fread (고속)" else "read.csv"))
lg(sprintf("병렬 코어: %d개  배치 크기: %d파일", N_CORES, BATCH_SIZE))

# =============================================================================
# ■ 테이블 1 & 2: 문항 1 응답 분포
# =============================================================================
resp_files <- list.files(
  "output/responses", pattern = "_response\\.csv$", full.names = TRUE
)

val_resp <- validate_files(resp_files, "응답(response)")

if (length(resp_files) == 0) {
  lg("응답 파일 없음 — 테이블 1·2 생성 생략")
} else {
  lg_section("테이블 1·2: 문항 1 응답 분포 처리")
  n_files <- length(resp_files)

  # ── 단일 파일 처리 함수 (문항 1만) ─────────────────────────────────────────
  process_resp_item1 <- function(path, read_fn, n_cat, iv1, iv2, iv3, iv4,
                                 decode_fn, parse_fn) {
    meta <- parse_fn(path)
    iv   <- decode_fn(meta$cond_code)

    df <- tryCatch(read_fn(path), error = function(e) NULL)
    if (is.null(df)) return(NULL)

    # 유령 응답자(theta=NA) 제외
    df <- df[!is.na(df$theta), , drop = FALSE]
    n_persons <- nrow(df)
    if (n_persons == 0L || !"item1" %in% names(df)) return(NULL)

    vals   <- as.integer(df[["item1"]])
    counts <- tabulate(vals + 1L, nbins = n_cat)   # 0→bin1 ~ 4→bin5
    props  <- counts / n_persons

    list(
      # Table 1 행 (character vector)
      row = c(meta$cond_code, meta$rep_id,
              iv$iv1_sf4, iv$iv2_b_interval, iv$iv3_b_mean, iv$iv4_theta_dist,
              n_persons,
              counts, round(props, 6)),
      # Table 2 accumulator 갱신용
      cond_code      = meta$cond_code,
      iv1_sf4        = iv$iv1_sf4,
      iv2_b_interval = iv$iv2_b_interval,
      iv3_b_mean     = iv$iv3_b_mean,
      iv4_theta_dist = iv$iv4_theta_dist,
      counts         = counts,
      props          = props,
      n_persons      = n_persons
    )
  }

  # Table 1: 헤더
  T1_PATH <- "output/analysis/resp_dist_rep_item1.csv"
  writeLines(paste(c("cond_code","rep_id","iv1_sf4","iv2_b_interval",
                     "iv3_b_mean","iv4_theta_dist","n_persons",
                     paste0("n", 0:4), paste0("prop", 0:4)),
                   collapse = ","), T1_PATH)

  # Table 2: accumulator (key = cond_code)
  acc2 <- new.env(hash = TRUE, parent = emptyenv())

  total_t1_rows <- 0L
  n_batches <- ceiling(n_files / BATCH_SIZE)

  for (b in seq_len(n_batches)) {
    idx_s <- (b - 1L) * BATCH_SIZE + 1L
    idx_e <- min(b * BATCH_SIZE, n_files)
    batch <- resp_files[idx_s:idx_e]

    if (N_CORES > 1L) {
      cl <- makeCluster(N_CORES)
      clusterExport(cl,
        c("read_fast","USE_FREAD","N_CAT","IV1_MAP","IV2_MAP","IV3_MAP","IV4_MAP",
          "decode_cond","parse_fname","process_resp_item1"),
        envir = environment())
      results <- parLapply(cl, batch, function(p)
        process_resp_item1(p, read_fast, N_CAT, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                           decode_cond, parse_fname))
      stopCluster(cl)
    } else {
      results <- lapply(batch, function(p)
        process_resp_item1(p, read_fast, N_CAT, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                           decode_cond, parse_fname))
    }

    results <- results[!sapply(results, is.null)]
    if (length(results) == 0L) { gc(verbose = FALSE); next }

    # Table 1 append
    t1_mat <- do.call(rbind, lapply(results, `[[`, "row"))
    write.table(t1_mat, T1_PATH, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    total_t1_rows <- total_t1_rows + nrow(t1_mat)

    # Table 2 accumulator 갱신
    for (r in results) {
      key <- r$cond_code
      if (exists(key, envir = acc2, inherits = FALSE)) {
        old <- get(key, envir = acc2)
        old$sum_n       <- old$sum_n       + r$counts
        old$sum_prop    <- old$sum_prop    + r$props
        old$sum_persons <- old$sum_persons + r$n_persons
        old$count       <- old$count       + 1L
        assign(key, old, envir = acc2)
      } else {
        assign(key, list(
          cond_code      = r$cond_code,
          iv1_sf4        = r$iv1_sf4,
          iv2_b_interval = r$iv2_b_interval,
          iv3_b_mean     = r$iv3_b_mean,
          iv4_theta_dist = r$iv4_theta_dist,
          sum_n = r$counts, sum_prop = r$props,
          sum_persons = r$n_persons, count = 1L
        ), envir = acc2)
      }
    }

    rm(results, t1_mat); gc(verbose = FALSE)
    cat(sprintf("\r  응답 진행: %d/%d 배치 (%d/%d 파일)...",
                b, n_batches, idx_e, n_files))
  }
  cat("\n")

  lg(sprintf("  저장 완료: %s (%d행)", T1_PATH, total_t1_rows))

  # Table 2: 평균 계산
  cond2_df <- do.call(rbind, lapply(ls(envir = acc2), function(key) {
    a <- get(key, envir = acc2)
    data.frame(
      cond_code      = a$cond_code,
      iv1_sf4        = a$iv1_sf4,
      iv2_b_interval = a$iv2_b_interval,
      iv3_b_mean     = a$iv3_b_mean,
      iv4_theta_dist = a$iv4_theta_dist,
      n_reps         = a$count,
      mean_n_persons = round(a$sum_persons / a$count, 2),
      mean_n0 = a$sum_n[1]/a$count, mean_n1 = a$sum_n[2]/a$count,
      mean_n2 = a$sum_n[3]/a$count, mean_n3 = a$sum_n[4]/a$count,
      mean_n4 = a$sum_n[5]/a$count,
      mean_prop0 = round(a$sum_prop[1]/a$count, 6),
      mean_prop1 = round(a$sum_prop[2]/a$count, 6),
      mean_prop2 = round(a$sum_prop[3]/a$count, 6),
      mean_prop3 = round(a$sum_prop[4]/a$count, 6),
      mean_prop4 = round(a$sum_prop[5]/a$count, 6),
      stringsAsFactors = FALSE
    )
  }))
  cond2_df <- cond2_df[order(cond2_df$cond_code), ]

  T2_PATH <- "output/analysis/resp_dist_cond_item1.csv"
  write.csv(cond2_df, T2_PATH, row.names = FALSE)
  lg(sprintf("  저장 완료: %s (%d행)", T2_PATH, nrow(cond2_df)))

  # 응시자 수 분포 요약 로그
  np_vals <- sapply(ls(envir = acc2), function(k) {
    a <- get(k, envir = acc2); a$sum_persons / a$count
  })
  lg(sprintf("  응시자 수(mean_n_persons) 분포: 최솟값=%.0f  최댓값=%.0f  평균=%.1f",
             min(np_vals), max(np_vals), mean(np_vals)))

  rm(acc2); gc(verbose = FALSE)
}

# =============================================================================
# ■ 테이블 3: 문항 1 경계모수 쌍 전치 정리 (est_params 파일 기반)
# =============================================================================
est_files <- list.files(
  "output/estimated_params", pattern = "_est_params\\.csv$", full.names = TRUE
)

val_est <- validate_files(est_files, "추정파라미터(est_params)")

if (length(est_files) == 0) {
  lg("추정 파라미터 파일 없음 — 테이블 3 생성 생략")
} else {
  lg_section("테이블 3: 문항 1 경계모수 쌍 전치 처리")
  n_est     <- length(est_files)
  n_bat_est <- ceiling(n_est / BATCH_SIZE)

  USE_FREAD3 <- requireNamespace("data.table", quietly = TRUE)
  read_est   <- if (USE_FREAD3)
    function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
  else
    function(p) read.csv(p, stringsAsFactors = FALSE)

  # ── 단일 est_params 파일 처리 ───────────────────────────────────────────────
  process_est_item1 <- function(path, read_fn, iv1, iv2, iv3, iv4,
                                decode_fn, parse_fn) {
    meta <- parse_fn(path)
    iv   <- decode_fn(meta$cond_code)

    df <- tryCatch(read_fn(path), error = function(e) NULL)
    if (is.null(df)) return(NULL)

    # 추정 실패 파일
    if ("success" %in% names(df) && !isTRUE(as.logical(df$success[1]))) {
      return(list(
        row = c(meta$cond_code, meta$rep_id,
                iv$iv1_sf4, iv$iv2_b_interval, iv$iv3_b_mean, iv$iv4_theta_dist,
                "FALSE", NA, NA, NA, NA, NA, NA, NA, NA),
        cond_code = meta$cond_code, iv = iv,
        conv = FALSE, tany = NA_integer_,
        t12 = NA_integer_, t23 = NA_integer_, t34 = NA_integer_
      ))
    }

    item1 <- df[df$item %in% c("Item1","item1",1L), ]
    if (nrow(item1) == 0L) item1 <- df[1L, ]

    b1   <- item1$b1[1]; b2 <- item1$b2[1]
    b3   <- item1$b3[1]; b4 <- item1$b4[1]
    conv <- isTRUE(as.logical(item1$converged[1]))

    if (conv && !anyNA(c(b1, b2, b3, b4))) {
      t12  <- as.integer(b1 >= b2); t23 <- as.integer(b2 >= b3)
      t34  <- as.integer(b3 >= b4); tany <- as.integer(t12 | t23 | t34)
    } else {
      t12 <- t23 <- t34 <- tany <- NA_integer_
    }

    list(
      row = c(meta$cond_code, meta$rep_id,
              iv$iv1_sf4, iv$iv2_b_interval, iv$iv3_b_mean, iv$iv4_theta_dist,
              as.character(conv), b1, b2, b3, b4, tany, t12, t23, t34),
      cond_code = meta$cond_code, iv = iv,
      conv = conv, tany = tany, t12 = t12, t23 = t23, t34 = t34
    )
  }

  # Table 3a 헤더
  T3_PATH <- "output/analysis/transposition_pairs_rep.csv"
  writeLines(paste(c("cond_code","rep_id","iv1_sf4","iv2_b_interval",
                     "iv3_b_mean","iv4_theta_dist","converged",
                     "b1","b2","b3","b4","trans_any","trans_12","trans_23","trans_34"),
                   collapse = ","), T3_PATH)

  # Table 3b accumulator (key = cond_code)
  acc3 <- new.env(hash = TRUE, parent = emptyenv())
  total_t3_rows <- 0L

  for (b in seq_len(n_bat_est)) {
    idx_s <- (b - 1L) * BATCH_SIZE + 1L
    idx_e <- min(b * BATCH_SIZE, n_est)
    batch <- est_files[idx_s:idx_e]

    if (N_CORES > 1L) {
      cl <- makeCluster(N_CORES)
      clusterExport(cl,
        c("read_est","USE_FREAD3","IV1_MAP","IV2_MAP","IV3_MAP","IV4_MAP",
          "decode_cond","parse_fname","process_est_item1"),
        envir = environment())
      results <- parLapply(cl, batch, function(p)
        process_est_item1(p, read_est, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                          decode_cond, parse_fname))
      stopCluster(cl)
    } else {
      results <- lapply(batch, function(p)
        process_est_item1(p, read_est, IV1_MAP, IV2_MAP, IV3_MAP, IV4_MAP,
                          decode_cond, parse_fname))
    }

    results <- results[!sapply(results, is.null)]
    if (length(results) == 0L) { gc(verbose = FALSE); next }

    # Table 3a append
    t3_mat <- do.call(rbind, lapply(results, `[[`, "row"))
    write.table(t3_mat, T3_PATH, append = TRUE, sep = ",",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    total_t3_rows <- total_t3_rows + nrow(t3_mat)

    # Table 3b accumulator 갱신
    for (r in results) {
      key <- r$cond_code
      cv  <- as.integer(isTRUE(r$conv))
      if (exists(key, envir = acc3, inherits = FALSE)) {
        old <- get(key, envir = acc3)
        old$n_reps <- old$n_reps + 1L
        old$n_conv <- old$n_conv + cv
        old$n_any  <- old$n_any  + ifelse(is.na(r$tany), 0L, r$tany)
        old$n_12   <- old$n_12   + ifelse(is.na(r$t12),  0L, r$t12)
        old$n_23   <- old$n_23   + ifelse(is.na(r$t23),  0L, r$t23)
        old$n_34   <- old$n_34   + ifelse(is.na(r$t34),  0L, r$t34)
        assign(key, old, envir = acc3)
      } else {
        assign(key, list(
          cond_code = r$cond_code, iv = r$iv,
          n_reps = 1L, n_conv = cv,
          n_any = ifelse(is.na(r$tany), 0L, r$tany),
          n_12  = ifelse(is.na(r$t12),  0L, r$t12),
          n_23  = ifelse(is.na(r$t23),  0L, r$t23),
          n_34  = ifelse(is.na(r$t34),  0L, r$t34)
        ), envir = acc3)
      }
    }

    rm(results, t3_mat); gc(verbose = FALSE)
    cat(sprintf("\r  추정 진행: %d/%d 배치 (%d/%d 파일)...",
                b, n_bat_est, idx_e, n_est))
  }
  cat("\n")
  lg(sprintf("  저장 완료: %s (%d행)", T3_PATH, total_t3_rows))

  # Table 3b 집계
  pairs_cond <- do.call(rbind, lapply(ls(envir = acc3), function(key) {
    a  <- get(key, envir = acc3)
    nc <- a$n_conv
    pct <- function(x) if (nc > 0L) round(100 * x / nc, 2) else NA_real_
    data.frame(
      cond_code      = a$cond_code,
      iv1_sf4        = a$iv$iv1_sf4,
      iv2_b_interval = a$iv$iv2_b_interval,
      iv3_b_mean     = a$iv$iv3_b_mean,
      iv4_theta_dist = a$iv$iv4_theta_dist,
      n_reps         = a$n_reps,
      n_converged    = nc,
      pct_converged  = round(100 * nc / a$n_reps, 2),
      n_trans_any    = a$n_any, pct_trans_any = pct(a$n_any),
      n_trans_12     = a$n_12,  pct_trans_12  = pct(a$n_12),
      n_trans_23     = a$n_23,  pct_trans_23  = pct(a$n_23),
      n_trans_34     = a$n_34,  pct_trans_34  = pct(a$n_34),
      stringsAsFactors = FALSE
    )
  }))
  pairs_cond <- pairs_cond[order(pairs_cond$cond_code), ]

  T3B_PATH <- "output/analysis/transposition_pairs_cond.csv"
  write.csv(pairs_cond, T3B_PATH, row.names = FALSE)
  lg(sprintf("  저장 완료: %s (%d행)", T3B_PATH, nrow(pairs_cond)))

  rm(acc3); gc(verbose = FALSE)
}

# =============================================================================
# 최종 요약 및 로그 저장
# =============================================================================
lg_section("최종 요약")

if (!is.null(val_resp)) {
  lg(sprintf("  [응답 파일]  발견 조건 %d개  누락 %d개  최대 반복 %d회",
             length(val_resp$reps_tbl), length(val_resp$missing_conds),
             val_resp$max_reps))
}
if (!is.null(val_est)) {
  lg(sprintf("  [추정 파일]  발견 조건 %d개  누락 %d개  최대 반복 %d회",
             length(val_est$reps_tbl), length(val_est$missing_conds),
             val_est$max_reps))
}
if (exists("total_t1_rows"))
  lg(sprintf("  테이블 1 (resp_dist_rep_item1)   : %d행", total_t1_rows))
if (exists("cond2_df"))
  lg(sprintf("  테이블 2 (resp_dist_cond_item1)  : %d행", nrow(cond2_df)))
if (exists("total_t3_rows"))
  lg(sprintf("  테이블 3a (transposition_pairs_rep) : %d행", total_t3_rows))
if (exists("pairs_cond"))
  lg(sprintf("  테이블 3b (transposition_pairs_cond): %d행", nrow(pairs_cond)))

lg(sprintf("\n실행 종료: %s", Sys.time()))
lg("================================================================")

flush_log()
cat(sprintf("\n로그 저장 완료: %s\n", LOG_PATH))
