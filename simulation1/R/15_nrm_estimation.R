# =============================================================================
# 15_nrm_estimation.R  (simulation1)
# 06_main.R 이 생성한 응답에 NRM(명명반응모형)을 적용하고,
# 문항 1의 경계모수 b1~b4 와 그 SE, b3·b4 차이 검정을 반복마다 계산
# (simulation2 의 11_nrm_estimation.R 과 같은 절차·같은 출력 열 + SE 열 추가)
#
# ── 1) NRM 추정 ───────────────────────────────────────────────────────────────
#   mirt itemtype = "nominal" (SE = TRUE), 20문항 전체를 함께 추정
#   P(X=k|θ) ∝ exp(a1·ak_k·θ + d_k),   k = 0..4
#   식별 제약(mirt 기본): ak0 = 0, ak4 = 4, d0 = 0  → 추정 모수: a1, ak1~ak3, d1~d4
#
# ── 2) 경계모수 (인접 두 범주의 응답 확률이 같아지는 θ) ──────────────────────
#   b_k = (d_{k-1} − d_k) / (a1 · (ak_k − ak_{k-1})),   k = 1..4
#   b3 = 범주 3과 4 사이 (범주 index 1부터), b4 = 범주 4와 5 사이
#
# ── 3) SE (델타 방법, mirt 공분산 행렬 사용) ─────────────────────────────────
#   se_a1, se_ak1~se_ak3, se_d1~se_d4 : NRM 추정 모수의 SE (공분산 행렬 대각)
#   se_b1~se_b4                      : 경계모수 SE  = sqrt(∇b_kᵀ V ∇b_k)
#   se_diff                          : SE(b4 − b3)  = sqrt(∇(b4−b3)ᵀ V ∇(b4−b3))
#   (기울기 ∇ 는 수치 미분, V 는 문항 1 자유모수의 공분산 행렬)
#
# ── 4) 양측 검정 ─────────────────────────────────────────────────────────────
#   H0: b4 = b3,  z = (b4 − b3) / se_diff,  |z| > Z_CRIT(1.96) 이면 기각, p = 2·P(Z > |z|)
#
# ── 5) 채점함수 전치 ──────────────────────────────────────────────────────────
#   sf_rev_k = (ak_k ≤ ak_{k-1}),  k = 1..4  (b_k 분모와 같은 위치), sf_transposed_any = 하나라도
#
# ── 척도 주의 (균등분포 조건) ────────────────────────────────────────────────
#   mirt 는 θ ~ N(0, 1) 로 척도를 고정하므로 추정값은 평균 0, 분산 1 척도에 있음
#   정규·편포 조건은 생성 분포도 평균 0, 분산 1 → 생성 척도와 같음
#   균등분포 U(−4, 4) 는 표준편차 8/√12 ≈ 2.31 → 참값과 비교할 때 b / 2.31 로 바꿔야 함
#   (13_analysis_limit_bias.R 의 b_true_std 와 같은 처리)
#
# 입력 : output/responses/*_response.csv  (조건 코드 4자리, 06_main.R 출력 그대로)
# 출력 : output/nrm/estimated_params/<응답파일명>_nrm.csv  (반복당 1행)
#        output/nrm/15_log.txt
#
# 이미 결과 파일이 있는 반복은 건너뜀 → 중단 후 다시 실행하면 이어서 진행
# 실행: simulation1 폴더를 작업 디렉토리로 두고 source("R/15_nrm_estimation.R")
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
COND_CODES_OVERRIDE <- NULL   # NULL = 응답 파일이 있는 모든 조건. 예: c("3112", "3114")
THETA_DIST_DIGITS   <- NULL   # 조건 코드 넷째 자리로 능력 분포 선택. NULL = 전부
                              # 1 = 정적편포, 2 = 정규, 3 = 부적편포, 4 = 균등
                              # 예: 정규분포만 → "2"
USE_PARALLEL <- TRUE
N_CORES      <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
BATCH_SIZE   <- 200L          # 배치당 파일 수 (진행 표시·메모리 관리용)
Z_CRIT       <- 1.96          # 양측 임계값 (α = .05)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(parallel)
  library(mirt)
})

RESP_DIR <- "output/responses"
NRM_DIR  <- "output/nrm/estimated_params"
dir.create(NRM_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- "output/nrm/15_log.txt"
lg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_PATH, append = TRUE)
}

# ── NRM 모수 → 경계모수 b1~b4 ─────────────────────────────────────────────────
nrm_boundaries <- function(p) {
  vapply(1:4, function(k) {
    (p[[paste0("d", k - 1)]] - p[[paste0("d", k)]]) /
      (p[["a1"]] * (p[[paste0("ak", k)]] - p[[paste0("ak", k - 1)]]))
  }, numeric(1))
}

# ── 반복 1개 처리: NRM 추정 → 경계모수·SE → 차이 검정 → CSV 저장 ────────────
fit_nrm_one <- function(resp_path, out_dir, z_crit) {
  base     <- sub("_response\\.csv$", "", basename(resp_path))
  out_path <- file.path(out_dir, paste0(base, "_nrm.csv"))

  par_names  <- c("a1", paste0("ak", 0:4), paste0("d", 0:4))
  free_names <- c("a1", paste0("ak", 1:3), paste0("d", 1:4))
  row <- c(
    list(cond_code = sub(".*_cond([0-9]+)_.*", "\\1", base),
         rep_id    = as.integer(sub(".*_rep([0-9]+)_.*", "\\1", base)),
         seed      = sub(".*_seed([0-9]+)$", "\\1", base),
         n_persons = NA_integer_, n_ghost = NA_integer_,
         success   = FALSE, converged = NA, se_ok = FALSE),
    setNames(as.list(rep(NA_real_, length(par_names))), par_names),
    setNames(as.list(rep(NA_real_, length(free_names))), paste0("se_", free_names)),
    list(b1 = NA_real_, b2 = NA_real_, b3 = NA_real_, b4 = NA_real_,
         se_b1 = NA_real_, se_b2 = NA_real_, se_b3 = NA_real_, se_b4 = NA_real_,
         diff_b4_b3 = NA_real_, se_diff = NA_real_,
         z = NA_real_, p_value = NA_real_, reject_h0 = NA,
         transposed_34 = NA, transposed_any = NA,
         sf_rev_1 = NA, sf_rev_2 = NA, sf_rev_3 = NA, sf_rev_4 = NA,
         sf_transposed_any = NA, error = "")
  )

  res <- tryCatch({
    df   <- read.csv(resp_path)
    resp <- as.data.frame(df[, grep("^item[0-9]+$", names(df))])
    row$n_persons <- nrow(resp)
    row$n_ghost   <- if ("theta" %in% names(df)) sum(is.na(df$theta)) else NA_integer_

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

    # 채점함수 전치: ak_k ≤ ak_{k-1}
    ak_rev <- diff(unlist(p[paste0("ak", 0:4)])) <= 0
    for (k in 1:4) row[[paste0("sf_rev_", k)]] <- ak_rev[k]
    row$sf_transposed_any <- any(ak_rev)

    # ── SE: 공분산 행렬 (mirt 추출 함수 사용 — stats::vcov 는 mirt S4 메서드로 연결되지 않음)
    #    error 메시지는 Windows 인코딩 깨짐 방지를 위해 영어로 기록
    se_msg <- ""
    V <- tryCatch(mirt::extract.mirt(mod, "vcov"), error = function(e) {
      se_msg <<- paste("vcov extract failed:", conditionMessage(e)); NULL
    })
    free_nm <- it1$name[it1$est]
    free_pn <- it1$parnum[it1$est]
    if (!is.null(V)) {
      rn  <- rownames(V)
      idx <- match(paste0(free_nm, ".", free_pn), rn)          # 예: "a1.1"
      if (anyNA(idx)) {
        idx <- match(free_pn, suppressWarnings(as.integer(sub(".*\\.", "", rn))))
      }
      if (anyNA(idx)) {
        se_msg <- paste0("vcov names not matched; first names: ",
                         paste(head(rn, 6), collapse = " "))
      } else {
        Vf <- V[idx, idx, drop = FALSE]

        # NRM 모수 SE
        dv     <- diag(Vf)
        se_par <- ifelse(is.finite(dv) & dv > 0, sqrt(pmax(dv, 0)), NA_real_)
        for (j in seq_along(free_nm))
          if (paste0("se_", free_nm[j]) %in% names(row)) row[[paste0("se_", free_nm[j])]] <- se_par[j]

        # 델타 방법: 함수 f(p) 의 분산 = ∇fᵀ V ∇f  (수치 미분)
        delta_var <- function(f) {
          grad <- vapply(free_nm, function(nm) {
            h  <- 1e-6 * max(1, abs(p[[nm]]))
            pu <- p; pu[[nm]] <- pu[[nm]] + h
            pd <- p; pd[[nm]] <- pd[[nm]] - h
            (f(pu) - f(pd)) / (2 * h)
          }, numeric(1))
          as.numeric(t(grad) %*% Vf %*% grad)
        }
        se_or_na <- function(v) if (is.finite(v) && v > 0) sqrt(v) else NA_real_

        for (k in 1:4)
          row[[paste0("se_b", k)]] <- se_or_na(delta_var(function(pp) nrm_boundaries(pp)[k]))

        var_d <- delta_var(function(pp) { bb <- nrm_boundaries(pp); bb[4] - bb[3] })
        if (is.finite(var_d) && var_d > 0) {
          row$se_ok     <- TRUE
          row$se_diff   <- sqrt(var_d)
          row$z         <- row$diff_b4_b3 / row$se_diff
          row$p_value   <- 2 * pnorm(abs(row$z), lower.tail = FALSE)
          row$reject_h0 <- abs(row$z) > z_crit
        } else {
          se_msg <- sprintf("invalid variance of b4-b3: %s (vcov has NA: %s)",
                            format(var_d), anyNA(Vf))
        }
      }
    } else if (se_msg == "") {
      se_msg <- "vcov is NULL (information matrix not computed)"
    }
    if (!row$se_ok) row$error <- se_msg
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
resp_files <- list.files(RESP_DIR, pattern = "_cond[0-9]{4}_rep[0-9]+_.*_response\\.csv$",
                         full.names = TRUE)
if (length(resp_files) == 0) stop("응답 파일이 없습니다: ", RESP_DIR)

fc <- sub(".*_cond([0-9]+)_.*", "\\1", basename(resp_files))
if (!is.null(COND_CODES_OVERRIDE)) {
  resp_files <- resp_files[fc %in% COND_CODES_OVERRIDE]; fc <- fc[fc %in% COND_CODES_OVERRIDE]
}
if (!is.null(THETA_DIST_DIGITS)) {
  resp_files <- resp_files[substr(fc, 4, 4) %in% THETA_DIST_DIGITS]
}
if (length(resp_files) == 0) stop("선택한 조건에 해당하는 응답 파일이 없습니다.")

done_base  <- sub("_nrm\\.csv$", "", list.files(NRM_DIR, pattern = "_nrm\\.csv$"))
todo_files <- resp_files[!sub("_response\\.csv$", "", basename(resp_files)) %in% done_base]

lg("================================================================")
lg(sprintf("NRM 추정 시작 (simulation1) — 응답 파일 %d개 (조건 %d개), 완료 %d개, 남은 파일 %d개",
           length(resp_files),
           length(unique(sub(".*_cond([0-9]+)_.*", "\\1", basename(resp_files)))),
           length(resp_files) - length(todo_files), length(todo_files)))
lg(sprintf("병렬: %s (%d 코어), Z_CRIT = %.2f", USE_PARALLEL, N_CORES, Z_CRIT))

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

lg("NRM 추정 완료 → 결과: output/nrm/estimated_params/")
