# =============================================================================
# 01_collect_data.R  (simulation1/R_graph)
# 그래프에 쓸 문항 1 추정 결과를 모아 한 번만 저장 (이후 그래프 스크립트는 이 파일만 읽음)
#
# 시뮬레이션1 (GPCM 추정): 정규분포 27개 조건
#   입력: output/estimated_params/*_cond[C]_rep[R]_seed[S]_est_params.csv  (C 넷째 자리 = 2)
#
# 시뮬레이션2 (NRM 추정): 시뮬레이션1과 같은 27개 조건만 (임시)
#   입력: ../simulation2/output/nrm/estimated_params/*_nrm.csv
#   경계모수 간격 1.5 는 시뮬레이션2에 없으므로 1.4 결과를 가져옴 (00_config.R 의 SIM2_SUBSTITUTE_IV2)
#     interval      = 그래프에 쓰는 수준 (1.5)
#     interval_true = 실제 생성에 쓴 수준 (1.4) → 참값 계산에 사용
#
# 출력: R_graph/output/data/
#   sim1_gpcm_item1.rds — 반복별 문항 1 GPCM 추정값 (a, b1~b4) 과 참값
#   sim2_nrm_item1.rds  — 반복별 문항 1 NRM 추정값 (a1, ak0~ak4, d0~d4, b1~b4) 과 참값
#   collect_log.txt     — 조건별 반복 수, 누락 조건
#
# 결과 파일이 새로 생기면 이 스크립트만 다시 실행
# =============================================================================

source("R_graph/00_functions.R", encoding = "UTF-8")

log_lines <- character(0)
lg <- function(...) { msg <- paste0(...); message(msg); log_lines <<- c(log_lines, msg) }

read_one <- if (requireNamespace("data.table", quietly = TRUE)) {
  function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE)
} else {
  function(p) read.csv(p, stringsAsFactors = FALSE)
}

# 같은 조건·반복 파일이 여러 개면 (다시 실행한 경우) 파일명 순서상 마지막(최신 날짜)만 사용
dedupe_files <- function(files) {
  files <- sort(files)
  key   <- sub("^[0-9]+_", "", sub("_seed[0-9]+.*$", "", basename(files)))
  dup   <- duplicated(key, fromLast = TRUE)
  if (any(dup)) lg(sprintf("  [주의] 같은 조건·반복 파일 %d개 중복 → 최신 파일만 사용", sum(dup)))
  files[!dup]
}

read_files <- function(files, reader) {
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    out[[i]] <- reader(files[i])
    if (i %% 2000 == 0) message(sprintf("    %d / %d", i, length(files)))
  }
  bind_rows(out)
}

check_counts <- function(df, label, expected_codes) {
  cnt <- df %>% count(cond_code, name = "n_rep")
  miss <- setdiff(expected_codes, cnt$cond_code)
  lg(sprintf("  %s: %d개 조건, 반복 수 %s", label, nrow(cnt),
             if (nrow(cnt)) paste(range(cnt$n_rep), collapse = " ~ ") else "0"))
  if (length(miss)) lg(sprintf("  [경고] %s 누락 조건: %s", label, paste(miss, collapse = ", ")))
}

grid <- expand.grid(mean = SIM1_IV3, interval = SIM1_IV2, sf4 = SIM1_IV1)[, 3:1]

# =============================================================================
# 1. 시뮬레이션1 (GPCM)
# =============================================================================
lg("[시뮬레이션1] GPCM 추정 결과")
grid$code1 <- with(grid, sim1_code(sf4, interval, mean))

est_dir <- file.path(SIM1_OUTPUT, "estimated_params")
files1  <- list.files(est_dir, pattern = "_est_params\\.csv$", full.names = TRUE)
code_f1 <- sub(".*_cond([0-9]+)_.*", "\\1", basename(files1))
files1  <- dedupe_files(files1[code_f1 %in% grid$code1])
if (length(files1) == 0) stop("시뮬레이션1 결과 파일이 없습니다: ", normalizePath(est_dir, mustWork = FALSE))
lg(sprintf("  파일 %d개 읽는 중...", length(files1)))

read_gpcm <- function(path) {
  code <- sub(".*_cond([0-9]+)_.*", "\\1", basename(path))
  rep  <- as.integer(sub(".*_rep([0-9]+)_.*", "\\1", basename(path)))
  df   <- tryCatch(read_one(path), error = function(e) NULL)
  na_row <- data.frame(cond_code = code, rep_id = rep, converged = FALSE,
                       a = NA_real_, b1 = NA_real_, b2 = NA_real_, b3 = NA_real_, b4 = NA_real_)
  if (is.null(df) || !"b1" %in% names(df)) return(na_row)
  it <- df[tolower(as.character(df$item)) %in% c("item1", "1"), , drop = FALSE]
  if (nrow(it) == 0) it <- df[1, , drop = FALSE]
  data.frame(cond_code = code, rep_id = rep,
             converged = as.logical(toupper(as.character(it$converged[1]))),
             a = it$a[1], b1 = it$b1[1], b2 = it$b2[1], b3 = it$b3[1], b4 = it$b4[1])
}

sim1 <- read_files(files1, read_gpcm) %>%
  left_join(grid %>% select(cond_code = code1, sf4, interval, mean), by = "cond_code") %>%
  mutate(interval_true = interval)
check_counts(sim1, "시뮬레이션1", grid$code1)

# =============================================================================
# 2. 시뮬레이션2 (NRM) — 시뮬레이션1과 같은 조건만
# =============================================================================
lg("[시뮬레이션2] NRM 추정 결과")
sub_iv2 <- function(x) {
  if (!length(SIM2_SUBSTITUTE_IV2)) return(x)
  s <- unname(SIM2_SUBSTITUTE_IV2[fmt_level(x)])
  ifelse(is.na(s), x, s)
}
grid <- grid %>%
  mutate(interval2 = sub_iv2(interval),
         i1 = match(sf4, SIM2_IV1), i2 = match(interval2, SIM2_IV2), i3 = match(mean, SIM2_IV3),
         code2 = sprintf("%02d%02d%02d%02d", i1, i2, i3, SIM2_IV4_INDEX))
if (anyNA(grid[, c("i1", "i2", "i3")])) {
  bad <- grid[!complete.cases(grid[, c("i1", "i2", "i3")]), c("sf4", "interval2", "mean")]
  stop("시뮬레이션2에 없는 수준이 있습니다 (SIM2_SUBSTITUTE_IV2 확인):\n",
       paste(capture.output(print(bad)), collapse = "\n"))
}
if (length(SIM2_SUBSTITUTE_IV2))
  lg(sprintf("  임시 대체: 경계모수 간격 %s",
             paste(names(SIM2_SUBSTITUTE_IV2), "←", SIM2_SUBSTITUTE_IV2, collapse = ", ")))

nrm_dir <- file.path(SIM2_OUTPUT, "nrm", "estimated_params")
files2  <- list.files(nrm_dir, pattern = "_nrm\\.csv$", full.names = TRUE)
code_f2 <- sub(".*_cond([0-9]+)_.*", "\\1", basename(files2))
files2  <- dedupe_files(files2[code_f2 %in% grid$code2])

if (length(files2) == 0) {
  n_all <- length(code_f2)
  if (!dir.exists(nrm_dir)) {
    lg(sprintf("  [경고] 시뮬레이션2 NRM 결과 폴더가 없습니다: %s", normalizePath(nrm_dir, mustWork = FALSE)))
  } else if (n_all == 0) {
    lg(sprintf("  [경고] 폴더는 있지만 *_nrm.csv 파일이 없습니다: %s", normalizePath(nrm_dir)))
  } else {
    lg(sprintf("  [경고] *_nrm.csv 파일 %d개가 있지만 가져올 27개 조건과 맞는 조건 코드가 없습니다", n_all))
    lg(sprintf("         폴더에 있는 조건 코드 (앞 10개): %s",
               paste(head(sort(unique(code_f2)), 10), collapse = ", ")))
    lg(sprintf("         찾으려던 조건 코드 (앞 10개): %s", paste(head(sort(grid$code2), 10), collapse = ", ")))
    lg("         → 00_config.R 의 SIM2_IV1·SIM2_IV2·SIM2_IV3 를 시뮬레이션2 R/00_setup.R 의 수준 목록과 같게 맞추세요")
  }
  lg("         NRM 그래프(04·05·07·08번 일부)는 건너뜀")
  sim2 <- NULL
} else {
  lg(sprintf("  파일 %d개 읽는 중...", length(files2)))
  nrm_cols <- c("a1", paste0("ak", 0:4), paste0("d", 0:4), paste0("b", 1:4))
  read_nrm <- function(path) {
    code <- sub(".*_cond([0-9]+)_.*", "\\1", basename(path))
    rep  <- as.integer(sub(".*_rep([0-9]+)_.*", "\\1", basename(path)))
    df   <- tryCatch(read_one(path), error = function(e) NULL)
    row  <- data.frame(cond_code2 = code, rep_id = rep, converged = FALSE)
    for (nm in nrm_cols) row[[nm]] <- NA_real_
    if (is.null(df) || nrow(df) == 0) return(row)
    row$converged <- as.logical(toupper(as.character(df$converged[1]))) %in% TRUE &
                     as.logical(toupper(as.character(df$success[1]))) %in% TRUE
    for (nm in intersect(nrm_cols, names(df))) row[[nm]] <- suppressWarnings(as.numeric(df[[nm]][1]))
    row
  }
  # 조건 코드가 실제로 그 조건인지 시뮬레이션2 true_params 로 확인
  #   (수준 목록이 다르면 다른 조건의 코드와 우연히 맞을 수 있음)
  tp_dir   <- file.path(SIM2_OUTPUT, "true_params")
  tp_files <- list.files(tp_dir, pattern = "_item_params\\.csv$", full.names = TRUE)
  if (length(tp_files) == 0) {
    lg(sprintf("  [주의] %s 가 없어 조건 코드를 확인하지 못함 — 아래 대응표를 직접 확인하세요", tp_dir))
    for (i in seq_len(nrow(grid)))
      lg(sprintf("         s4 = %.2f, 간격 = %.2f, 난이도 = %.1f  →  cond%s",
                 grid$sf4[i], grid$interval2[i], grid$mean[i], grid$code2[i]))
  } else {
    tp_code <- sub(".*_cond([0-9]+)_.*", "\\1", basename(tp_files))
    bad <- character(0)
    for (i in seq_len(nrow(grid))) {
      f <- tp_files[tp_code == grid$code2[i]]
      if (length(f) == 0) next
      it <- read_one(f[1])[1, ]
      exp_b <- grid$mean[i] + grid$interval2[i] * B_OFFSET
      ok <- abs(it$sf4 - grid$sf4[i]) < 1e-6 &&
            all(abs(unlist(it[paste0("b", 1:4)]) - exp_b) < 1e-6)
      if (!ok) bad <- c(bad, sprintf("cond%s: 기대 s4 = %g, b = (%s) / 파일 s4 = %g, b = (%s)",
                                     grid$code2[i], grid$sf4[i], paste(exp_b, collapse = ", "),
                                     it$sf4, paste(round(unlist(it[paste0("b", 1:4)]), 3), collapse = ", ")))
    }
    if (length(bad))
      stop("시뮬레이션2 조건 코드와 실제 조건이 다릅니다. 00_config.R 의 SIM2_IV1·SIM2_IV2·SIM2_IV3 를\n",
           "시뮬레이션2 R/00_setup.R 의 수준 목록과 같게(순서까지) 맞추세요.\n  ",
           paste(head(bad, 5), collapse = "\n  "))
    lg("  조건 코드 확인 완료 (true_params 의 채점함수·경계모수와 일치)")
  }

  sim2 <- read_files(files2, read_nrm) %>%
    left_join(grid %>% select(cond_code2 = code2, cond_code = code1, sf4, interval,
                              interval_true = interval2, mean), by = "cond_code2")
  check_counts(sim2, "시뮬레이션2", grid$code1)
}

# =============================================================================
# 3. 참값 붙여 저장
# =============================================================================
add_true <- function(df) {
  for (k in 1:4) df[[paste0("b", k, "_true")]] <- df$mean + df$interval_true * B_OFFSET[k]
  df
}

sim1 <- add_true(sim1) %>% arrange(cond_code, rep_id)
saveRDS(sim1, file.path(DATA_DIR, "sim1_gpcm_item1.rds"))
lg(sprintf("저장: %s (%d행)", file.path(DATA_DIR, "sim1_gpcm_item1.rds"), nrow(sim1)))

if (!is.null(sim2)) {
  sim2 <- add_true(sim2) %>% arrange(cond_code, rep_id)
  saveRDS(sim2, file.path(DATA_DIR, "sim2_nrm_item1.rds"))
  lg(sprintf("저장: %s (%d행)", file.path(DATA_DIR, "sim2_nrm_item1.rds"), nrow(sim2)))
} else if (file.exists(file.path(DATA_DIR, "sim2_nrm_item1.rds"))) {
  file.remove(file.path(DATA_DIR, "sim2_nrm_item1.rds"))
}

writeLines(log_lines, file.path(DATA_DIR, "collect_log.txt"), useBytes = TRUE)
