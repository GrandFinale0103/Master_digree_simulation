# =============================================================================
# run_reanalysis.R  (simulation2)
# 전치 집계 수정 이후 다시 돌려야 하는 분석 스크립트를 한 번에 실행
#
# 사용법: 다른 스크립트와 같이 작업 디렉토리를 simulation2/ 로 둔 상태에서
#           source("R/run_reanalysis.R")
#         또는 터미널에서 simulation2/ 로 이동한 뒤
#           Rscript R/run_reanalysis.R
#
# 실행 순서 (괄호: 먼저 성공해야 하는 단계)
#   GPCM : 08_analysis_transposition.R → 10_plot.R (08), 13_plot_transposition_pairs.R (08)
#   NRM  : 12_analysis_nrm_test.R → 14_plot_nrm_transposition_pairs.R (12),
#          15_diagnose_nrm_extremes.R (12), 16_plot_nrm_excluded.R (12)
#   NRM 추정(11), 응답 분포(07), RMSE/Bias(09)는 다시 돌리지 않음 (기존 결과 사용)
#
# - 각 스크립트는 별도의 R 프로세스(Rscript)로 실행 → 스크립트 사이 변수 충돌 없음
# - 한 스크립트가 실패하면 그 결과를 쓰는 뒤 스크립트만 건너뛰고 나머지는 계속 실행
# - 스크립트별 출력: output/reanalysis_logs/<스크립트>.log
#   전체 요약     : output/reanalysis_logs/summary.txt
# =============================================================================

if (!file.exists("R/12_analysis_nrm_test.R"))
  stop("작업 디렉토리를 simulation2/ 로 맞춘 뒤 실행하세요. 현재: ", getwd())

LOG_DIR <- "output/reanalysis_logs"
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
RSCRIPT <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")

# 실행 순서와 의존 관계 (needs: 먼저 성공해야 하는 단계)
STEPS <- list(
  list(script = "08_analysis_transposition.R",       needs = character(0)),
  list(script = "10_plot.R",                         needs = "08_analysis_transposition.R"),
  list(script = "13_plot_transposition_pairs.R",     needs = "08_analysis_transposition.R"),
  list(script = "12_analysis_nrm_test.R",            needs = character(0)),
  list(script = "14_plot_nrm_transposition_pairs.R", needs = "12_analysis_nrm_test.R"),
  list(script = "15_diagnose_nrm_extremes.R",        needs = "12_analysis_nrm_test.R"),
  list(script = "16_plot_nrm_excluded.R",            needs = "12_analysis_nrm_test.R")
)

status <- setNames(rep(NA_character_, length(STEPS)), sapply(STEPS, `[[`, "script"))
summary_lines <- character(0)
note <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(msg)
  summary_lines <<- c(summary_lines, msg)
  writeLines(summary_lines, file.path(LOG_DIR, "summary.txt"))
}

note(sprintf("재분석 시작 — %d단계, 로그 폴더: %s", length(STEPS), normalizePath(LOG_DIR)))
total_start <- Sys.time()

for (st in STEPS) {
  failed_deps <- st$needs[status[st$needs] != "OK"]
  if (length(failed_deps) > 0) {
    status[st$script] <- "SKIP"
    note(sprintf("건너뜀  R/%s  (앞 단계 실패: %s)", st$script, paste(failed_deps, collapse = ", ")))
    next
  }

  log_path <- file.path(LOG_DIR, sub("\\.R$", ".log", st$script))
  note(sprintf("실행    R/%s", st$script))
  t0 <- Sys.time()

  code <- tryCatch(
    system2(RSCRIPT, c("--no-save", "--no-restore", shQuote(file.path("R", st$script))),
            stdout = log_path, stderr = log_path),
    error = function(e) { cat(conditionMessage(e), file = log_path, append = TRUE); -1L })

  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (identical(as.integer(code), 0L)) {
    status[st$script] <- "OK"
    note(sprintf("완료    R/%s (%.1f분)", st$script, mins))
  } else {
    status[st$script] <- "FAIL"
    tail_log <- tryCatch(tail(readLines(log_path, warn = FALSE), 8), error = function(e) character(0))
    note(sprintf("실패    R/%s (%.1f분, 종료 코드 %s) — 로그: %s", st$script, mins, code, log_path))
    for (l in tail_log) note("          | ", l)
  }
}

note(sprintf("전체 종료 — %.1f분 | 성공 %d, 실패 %d, 건너뜀 %d",
             as.numeric(difftime(Sys.time(), total_start, units = "mins")),
             sum(status == "OK"), sum(status == "FAIL"), sum(status == "SKIP")))
