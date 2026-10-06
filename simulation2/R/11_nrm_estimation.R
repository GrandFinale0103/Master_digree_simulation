# =============================================================================
# 11_nrm_estimation.R  (simulation2)
# 이미 생성된 응답 데이터에 NRM(명명반응모형)을 적용하고,
# 문항 1의 경계모수 b3, b4 차이에 대해 단방향 검정(Approach 2)을 반복마다 수행
#
# ── 1) NRM 추정 ───────────────────────────────────────────────────────────────
#   mirt itemtype = "nominal" (SE = TRUE)
#   P(X=k|θ) ∝ exp(a1·ak_k·θ + d_k),   k = 0..4
#   식별 제약(mirt 기본): ak0 = 0, ak4 = 4, d0 = 0
#
# ── 2) 경계모수 계산 (인접 두 범주의 응답 확률이 같아지는 θ) ──────────────────
#   a1·ak_{k-1}·θ + d_{k-1} = a1·ak_k·θ + d_k
#   ⇒ b_k = (d_{k-1} − d_k) / (a1 · (ak_k − ak_{k-1})),   k = 1..4
#   b3 = 범주 3과 4 사이 (범주 index 1부터), b4 = 범주 4와 5 사이
#
# ── 3) 단방향 검정 (Approach 2: 서열화된 경계모수 증거 탐색) ──────────────────
#   H0: b4 ≤ b3      HA: b4 > b3
#   z = (b4 − b3) / SE(b4 − b3),  SE는 델타 방법 (mirt 공분산 행렬 사용)
#   z > Z_CRIT(1.65) 이면 H0 기각 → "서열화됨(ordered)" 증거
#   p = P(Z > z)  (상단 단측)
#
# ── 4) 채점함수 전치 ──────────────────────────────────────────────────────────
#   추정된 채점함수 ak0..ak4 가 증가하지 않으면 전치로 표시
#   sf_rev_k = (ak_k ≤ ak_{k-1}),  k = 1..4  (b_k 계산의 분모와 같은 위치)
#   sf_transposed_any = 하나라도 전치
#
# 입력 : output/responses/*_response.csv  (06_main.R 이 생성한 응답 그대로 사용)
# 출력 : output/nrm/estimated_params/<응답파일명>_nrm.csv  (반복당 1행)
#        output/nrm/11_log.txt
#
# 이미 결과 파일이 있는 반복은 건너뜀 → 중단 후 다시 실행하면 이어서 진행
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
COND_CODES_OVERRIDE <- NULL   # NULL = 응답 파일이 있는 모든 조건
                              # 예: c("01010101", "07100401")
USE_PARALLEL <- TRUE
N_CORES      <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
BATCH_SIZE   <- 200L          # 배치당 파일 수 (진행 표시·메모리 관리용)
Z_CRIT       <- 1.65          # 단측 임계값
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(parallel)
  library(mirt)
})

RESP_DIR <- "output/responses"
NRM_DIR  <- "output/nrm/estimated_params"
dir.create(NRM_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- "output/nrm/11_log.txt"
lg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_PATH, append = TRUE)
}

# ── NRM 모수 → 경계모수 b1~b4 ─────────────────────────────────────────────────
nrm_boundaries <- function(p) {
  vapply(1:4, function(k) {
    (p[[paste0("d", k - 1)]] - p[[paste0("d", k)]]) /
      (p[[paste0("a1")]] * (p[[paste0("ak", k)]] - p[[paste0("ak", k - 1)]]))
  }, numeric(1))
}

# ── 반복 1개 처리: NRM 추정 → 경계모수 → 차이 검정 → CSV 저장 ────────────────
fit_nrm_one <- function(resp_path, out_dir, z_crit) {
  base     <- sub("_response\\.csv$", "", basename(resp_path))
  out_path <- file.path(out_dir, paste0(base, "_nrm.csv"))

  par_names <- c("a1", paste0("ak", 0:4), paste0("d", 0:4))
  row <- c(
    list(cond_code = sub(".*_cond([0-9]+)_.*", "\\1", base),
         rep_id    = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", base)),
         seed      = sub(".*_seed([0-9]+)$", "\\1", base),
         success   = FALSE, converged = NA, se_ok = FALSE),
    setNames(as.list(rep(NA_real_, length(par_names))), par_names),
    list(b1 = NA_real_, b2 = NA_real_, b3 = NA_real_, b4 = NA_real_,
         diff_b4_b3 = NA_real_, se_diff = NA_real_,
         z = NA_real_, p_value = NA_real_, reject_h0 = NA,
         transposed_34 = NA, transposed_any = NA,
         sf_rev_1 = NA, sf_rev_2 = NA, sf_rev_3 = NA, sf_rev_4 = NA,
         sf_transposed_any = NA, error = "")
  )

  res <- tryCatch({
    df   <- read.csv(resp_path)
    resp <- as.data.frame(df[, grep("^item[0-9]+$", names(df))])

    mod <- mirt::mirt(resp, model = 1, itemtype = "nominal", SE = TRUE,
                      verbose = FALSE, technical = list(NCYCLES = 3000))

    row$success   <- TRUE
    row$converged <- mirt::extract.mirt(mod, "converged")

    # 문항 1의 NRM 모수
    pv  <- mirt::mod2values(mod)
    it1 <- pv[pv$item == colnames(resp)[1], ]
    p   <- setNames(it1$value, it1$name)
    for (nm in par_names) row[[nm]] <- p[[nm]]

    b <- nrm_boundaries(p)
    row$b1 <- b[1]; row$b2 <- b[2]; row$b3 <- b[3]; row$b4 <- b[4]
    row$diff_b4_b3     <- b[4] - b[3]
    row$transposed_34  <- b[3] >= b[4]
    row$transposed_any <- any(diff(b) <= 0)

    # 채점함수 전치: ak_k ≤ ak_{k-1} (sf_rev_k = 범주 k와 k+1 사이, 1부터)
    ak_rev <- diff(unlist(p[paste0("ak", 0:4)])) <= 0
    for (k in 1:4) row[[paste0("sf_rev_", k)]] <- ak_rev[k]
    row$sf_transposed_any <- any(ak_rev)

    # 델타 방법: SE(b4 − b3)
    V <- tryCatch(stats::vcov(mod), error = function(e) NULL)
    free_nm <- it1$name[it1$est]
    free_pn <- it1$parnum[it1$est]
    if (!is.null(V)) {
      v_pn <- as.integer(sub(".*\\.", "", rownames(V)))
      idx  <- match(free_pn, v_pn)
      if (!anyNA(idx)) {
        Vf <- V[idx, idx, drop = FALSE]
        g  <- function(pp) { bb <- nrm_boundaries(pp); bb[4] - bb[3] }
        grad <- vapply(free_nm, function(nm) {
          h  <- 1e-6 * max(1, abs(p[[nm]]))
          pu <- p; pu[[nm]] <- pu[[nm]] + h
          pd <- p; pd[[nm]] <- pd[[nm]] - h
          (g(pu) - g(pd)) / (2 * h)
        }, numeric(1))
        var_d <- as.numeric(t(grad) %*% Vf %*% grad)
        if (is.finite(var_d) && var_d > 0) {
          row$se_ok     <- TRUE
          row$se_diff   <- sqrt(var_d)
          row$z         <- row$diff_b4_b3 / row$se_diff
          row$p_value   <- pnorm(row$z, lower.tail = FALSE)
          row$reject_h0 <- row$z > z_crit
        }
      }
    }
    if (!row$se_ok) row$error <- "SE 계산 불가 (정보행렬 역행렬 실패 등)"
    row
  }, error = function(e) {
    row$error <- conditionMessage(e)
    row
  })

  write.csv(as.data.frame(res, stringsAsFactors = FALSE), out_path, row.names = FALSE)
  res$success
}

# =============================================================================
# 실행
# =============================================================================
resp_files <- list.files(RESP_DIR, pattern = "_response\\.csv$", full.names = TRUE)
if (length(resp_files) == 0) stop("응답 파일이 없습니다: ", RESP_DIR)

if (!is.null(COND_CODES_OVERRIDE)) {
  fc <- sub(".*_cond([0-9]+)_.*", "\\1", basename(resp_files))
  resp_files <- resp_files[fc %in% COND_CODES_OVERRIDE]
}

done_base  <- sub("_nrm\\.csv$", "",
                  list.files(NRM_DIR, pattern = "_nrm\\.csv$"))
todo_files <- resp_files[!sub("_response\\.csv$", "", basename(resp_files)) %in% done_base]

lg("================================================================")
lg(sprintf("NRM 추정 시작 — 응답 파일 %d개, 완료 %d개, 남은 파일 %d개",
           length(resp_files), length(resp_files) - length(todo_files),
           length(todo_files)))
lg(sprintf("병렬: %s (%d 코어), Z_CRIT = %.2f",
           USE_PARALLEL, N_CORES, Z_CRIT))

if (length(todo_files) > 0) {
  start_time <- Sys.time()
  n_batches  <- ceiling(length(todo_files) / BATCH_SIZE)

  run_par <- USE_PARALLEL && N_CORES > 1
  if (run_par) {
    cl <- makeCluster(N_CORES)
    wd <- getwd()
    clusterExport(cl, c("nrm_boundaries", "fit_nrm_one", "wd"), envir = .GlobalEnv)
    clusterEvalQ(cl, { setwd(wd); suppressPackageStartupMessages(library(mirt)) })
  }

  n_fail <- 0L
  tryCatch({
    for (b in seq_len(n_batches)) {
      batch <- todo_files[((b - 1) * BATCH_SIZE + 1):min(b * BATCH_SIZE, length(todo_files))]
      ok <- if (run_par) {
        unlist(parLapply(cl, batch, fit_nrm_one, out_dir = NRM_DIR, z_crit = Z_CRIT))
      } else {
        vapply(batch, fit_nrm_one, logical(1), out_dir = NRM_DIR, z_crit = Z_CRIT)
      }
      n_fail <- n_fail + sum(!ok)
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      lg(sprintf("배치 %d/%d 완료 (누적 추정 실패 %d) — 경과 %.1f분",
                 b, n_batches, n_fail, elapsed))
    }
  }, finally = if (run_par) stopCluster(cl))
}

lg("NRM 추정 완료 → 다음 단계: source(\"R/12_analysis_nrm_test.R\")")
