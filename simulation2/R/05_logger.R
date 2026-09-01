# =============================================================================
# 05_logger.R
# 로그 기록 유틸리티
# =============================================================================

# ── 로그 항목 형식 ────────────────────────────────────────────────────────────
# [YYYY-MM-DD HH:MM:SS] rep=NNNN func=FUNC_NAME            status=OK   msg=...

# ── 로그 쓰기 ─────────────────────────────────────────────────────────────────
#
# parallel = FALSE : 로그 파일에 즉시 append
# parallel = TRUE  : output/logs/temp/ 에 rep별 임시 파일로 저장
#
write_log <- function(log_file, rep_id, func_name, status, msg = "", parallel = FALSE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- sprintf(
    "[%s] rep=%04d func=%-28s status=%-5s msg=%s\n",
    timestamp, rep_id, func_name, status, msg
  )

  if (parallel) {
    # 임시 파일명: gen_log_cond1111_rep0001.tmp 형식
    base    <- tools::file_path_sans_ext(basename(log_file))
    tmp_dir <- file.path(dirname(log_file), "temp")
    dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)
    tmp_file <- file.path(tmp_dir,
                          sprintf("%s_rep%04d.tmp", base, rep_id))
    cat(entry, file = tmp_file, append = TRUE)
  } else {
    cat(entry, file = log_file, append = TRUE)
  }
}

# ── 임시 로그 병합 (병렬 실행 완료 후 호출) ──────────────────────────────────
#
# output/logs/temp/ 의 .tmp 파일들을 rep 번호 순으로 정렬 후
# 대상 log_file 에 append, 완료 후 .tmp 파일 삭제
#
merge_temp_logs <- function(log_file) {
  base    <- tools::file_path_sans_ext(basename(log_file))
  tmp_dir <- file.path(dirname(log_file), "temp")

  tmp_files <- list.files(
    tmp_dir,
    pattern    = paste0("^", base, "_rep\\d+\\.tmp$"),
    full.names = TRUE
  )

  if (length(tmp_files) == 0) {
    cat("병합할 임시 로그 파일 없음:", log_file, "\n")
    return(invisible(NULL))
  }

  # rep 번호 기준 정렬
  tmp_files <- tmp_files[order(tmp_files)]

  for (f in tmp_files) {
    lines <- readLines(f, warn = FALSE)
    cat(paste(lines, collapse = "\n"), "\n", file = log_file, append = TRUE)
    file.remove(f)
  }

  cat(sprintf("로그 병합 완료: %d개 항목 → %s\n", length(tmp_files), log_file))
}

# ── 조건별 로그 파일 경로 반환 ───────────────────────────────────────────────
get_log_paths <- function(cond_code) {
  list(
    gen = file.path("output", "logs", sprintf("gen_log_cond%s.txt", cond_code)),
    est = file.path("output", "logs", sprintf("est_log_cond%s.txt", cond_code))
  )
}
