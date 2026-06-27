# =============================================================================
# 06_main.R
# 메인 실행 스크립트 — 여러 조건을 순서대로 실행
#
# 실행 전 준비사항:
#   1. 최초 1회: source("R/01_seed_generator.R") 으로 시드 테이블 생성
#   2. 아래 [사용자 설정 영역]의 COND_CODES 와 병렬 옵션 확인 후 실행
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                       사용자 설정 영역                                    ║
# ╠══════════════════════════════════════════════════════════════════════════╣
COND_CODES   <- c("1111", "1112", "1113", "1114")  # 실행할 조건 코드 목록
                         #   자리 1: 채점함수 sf[4] (1=3, 2=3.33, 3=3.66)
                         #   자리 2: 경계모수 간격  (1=2, 2=1.5,  3=1)
                         #   자리 3: 문항 심각도    (1=0, 2=1.5,  3=3)
                         #   자리 4: 능력모수 분포  (1=정적편포, 2=정규,
                         #                           3=부적편포, 4=균등)

USE_PARALLEL <- FALSE    # TRUE: 병렬 실행 / FALSE: 순차 실행
N_CORES      <- 4        # 병렬 시 사용 코어 수 (본인 환경에 맞게 조정)
N_REP        <- 1000     # 총 반복 횟수
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── 초기화 ───────────────────────────────────────────────────────────────────
source("R/00_setup.R")
source("R/01_seed_generator.R")   # seed_table.rds 없으면 자동 생성
source("R/02_item_params.R")
source("R/03_response.R")
source("R/04_estimation.R")
source("R/05_logger.R")

# ── 입력 검증 ────────────────────────────────────────────────────────────────
invisible(lapply(COND_CODES, parse_cond_code))  # 유효하지 않은 코드 있으면 에러

# ── 시드 테이블 로드 ─────────────────────────────────────────────────────────
seed_table <- readRDS("output/seed_table.rds")

# ── 단일 반복 실행 함수 ──────────────────────────────────────────────────────
run_one_rep <- function(rep_id, cond_code, cond_params, cond_seeds,
                        log_paths, progress_file, parallel) {

  seed <- cond_seeds$seed[cond_seeds$rep_id == rep_id]

  # ── 데이터 생성 단계 ────────────────────────────────────────────────────
  write_log(log_paths$gen, rep_id, "generate_item_params", "START",
            parallel = parallel)

  item_params <- tryCatch({
    ip <- generate_item_params(cond_params, rep_id, seed, cond_code)
    write_log(log_paths$gen, rep_id, "generate_item_params", "OK",
              parallel = parallel)
    ip
  }, error = function(e) {
    write_log(log_paths$gen, rep_id, "generate_item_params", "FAIL",
              msg = conditionMessage(e), parallel = parallel)
    NULL
  })

  if (is.null(item_params)) return(invisible(NULL))

  write_log(log_paths$gen, rep_id, "generate_response", "START",
            parallel = parallel)

  resp_data <- tryCatch({
    rd <- generate_response(item_params, cond_params, rep_id, seed, cond_code)
    write_log(log_paths$gen, rep_id, "generate_response", "OK",
              parallel = parallel)
    rd
  }, error = function(e) {
    write_log(log_paths$gen, rep_id, "generate_response", "FAIL",
              msg = conditionMessage(e), parallel = parallel)
    NULL
  })

  if (is.null(resp_data)) return(invisible(NULL))

  # ── 모수 추정 단계 ──────────────────────────────────────────────────────
  write_log(log_paths$est, rep_id, "estimate_params", "START",
            parallel = parallel)

  est_result <- tryCatch({
    er <- estimate_params(resp_data$response, rep_id, seed, cond_code)
    status <- if (isTRUE(er$success)) "OK" else "FAIL"
    msg    <- if (isTRUE(er$success)) {
      sprintf("converged=%s", er$convergence)
    } else {
      er$error
    }
    write_log(log_paths$est, rep_id, "estimate_params", status,
              msg = msg, parallel = parallel)
    er
  }, error = function(e) {
    write_log(log_paths$est, rep_id, "estimate_params", "FAIL",
              msg = conditionMessage(e), parallel = parallel)
    NULL
  })

  # ── 진행 기록 (순차 모드만: 병렬은 완료 후 일괄 기록) ──────────────────
  if (!parallel) {
    cat(rep_id, file = progress_file, append = TRUE, sep = "\n")
  }

  invisible(est_result)
}

# ── 조건별 실행 함수 ──────────────────────────────────────────────────────────
run_one_cond <- function(COND_CODE) {

  cond_params   <- parse_cond_code(COND_CODE)
  cond_seeds    <- seed_table[seed_table$cond_code == COND_CODE, ]
  progress_file <- file.path("output", "progress",
                             sprintf("cond%s_completed.txt", COND_CODE))
  log_paths     <- get_log_paths(COND_CODE)

  # 진행 상황 확인 (중단 후 재개 지원)
  completed_reps <- integer(0)
  if (file.exists(progress_file)) {
    completed_reps <- as.integer(readLines(progress_file, warn = FALSE))
    completed_reps <- completed_reps[!is.na(completed_reps)]
  }
  remaining_reps <- setdiff(seq_len(N_REP), completed_reps)

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat(sprintf("  조건 코드     : %s\n", COND_CODE))
  cat(sprintf("  채점함수 간격 : %.2f\n", cond_params$iv1_score_gap))
  cat(sprintf("  경계모수 간격 : %.2f\n", cond_params$iv2_b_interval))
  cat(sprintf("  문항 심각도   : %.2f\n", cond_params$iv3_b_mean))
  cat(sprintf("  능력모수 분포 : %s\n",   cond_params$iv4_theta_dist))
  cat(sprintf("  반복 횟수     : %d\n",   N_REP))
  cat(sprintf("  병렬 실행     : %s",     ifelse(USE_PARALLEL, "예", "아니오")))
  if (USE_PARALLEL) cat(sprintf(" (%d 코어)", N_CORES))
  cat("\n")
  cat(sprintf("  완료된 반복   : %d / %d\n", length(completed_reps), N_REP))
  cat(sprintf("  실행할 반복   : %d개\n", length(remaining_reps)))
  cat(strrep("=", 60), "\n\n", sep = "")

  if (length(remaining_reps) == 0) {
    cat("  → 이미 완료된 조건입니다. 건너뜁니다.\n\n")
    return(invisible(NULL))
  }

  start_time <- Sys.time()

  if (USE_PARALLEL) {

    cl <- makeCluster(N_CORES)
    registerDoParallel(cl)

    clusterExport(cl, varlist = c(
      "cond_params", "COND_CODE", "cond_seeds", "log_paths", "progress_file",
      "N_ITEMS", "N_CAT", "N_PERSONS", "DISCRIM",
      "IV1_LEVELS", "IV2_LEVELS", "IV3_LEVELS", "IV4_LEVELS",
      "RSN_POS", "RSN_NEG",
      "parse_cond_code", "generate_item_params",
      "gpcm_sf_prob", "sample_theta", "find_empty_cats", "generate_response",
      "estimate_params",
      "write_log", "get_log_paths"
    ), envir = environment())
    clusterEvalQ(cl, { library(mirt); library(sn) })

    completed_in_run <- foreach(
      rep_id = remaining_reps,
      .combine  = c,
      .packages = c("mirt", "sn")
    ) %dopar% {
      run_one_rep(rep_id, COND_CODE, cond_params, cond_seeds,
                  log_paths, progress_file, parallel = TRUE)
      rep_id
    }

    stopCluster(cl)

    cat(completed_in_run, file = progress_file, append = TRUE, sep = "\n")
    merge_temp_logs(log_paths$gen)
    merge_temp_logs(log_paths$est)

  } else {

    for (rep_id in remaining_reps) {
      if (rep_id %% 100 == 0 || rep_id == remaining_reps[1]) {
        cat(sprintf("[%s] 반복 %d / %d 진행 중...\n",
                    format(Sys.time(), "%H:%M:%S"), rep_id, N_REP))
      }
      run_one_rep(rep_id, COND_CODE, cond_params, cond_seeds,
                  log_paths, progress_file, parallel = FALSE)
    }
  }

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("\n조건 %s 완료 — 소요 시간: %.1f 분\n", COND_CODE, elapsed))
}

# ── 전체 실행 ────────────────────────────────────────────────────────────────
total_start <- Sys.time()
cat(sprintf("\n총 %d개 조건 실행 시작: %s\n",
            length(COND_CODES), paste(COND_CODES, collapse = ", ")))

for (COND_CODE in COND_CODES) {
  run_one_cond(COND_CODE)
}

total_elapsed <- difftime(Sys.time(), total_start, units = "mins")
cat(sprintf("\n전체 완료 — 총 소요 시간: %.1f 분\n", total_elapsed))
