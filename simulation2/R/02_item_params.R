# =============================================================================
# 02_item_params.R  (simulation2)
# 조건 코드 파싱 및 문항 모수 생성 함수
#
# 조건 코드 형식: 8자리 2자리 고정 (예: "01030201")
#   앞 2자리: IV1 인덱스 (01~07)
#   중 2자리: IV2 인덱스 (01~10)
#   중 2자리: IV3 인덱스 (01~04)
#   끝 2자리: IV4 인덱스 (01만 유효)
# =============================================================================

# ── 조건 코드 파싱 ────────────────────────────────────────────────────────────
parse_cond_code <- function(code) {
  if (!grepl("^[0-9]{8}$", code)) {
    stop(sprintf(
      "유효하지 않은 조건 코드: '%s' — 8자리 숫자여야 합니다 (예: '01030201')",
      code
    ))
  }
  d1 <- as.integer(substr(code, 1, 2))
  d2 <- as.integer(substr(code, 3, 4))
  d3 <- as.integer(substr(code, 5, 6))
  d4 <- as.integer(substr(code, 7, 8))

  n1 <- length(IV1_LEVELS); n2 <- length(IV2_LEVELS)
  n3 <- length(IV3_LEVELS); n4 <- length(IV4_LEVELS)

  if (d1 < 1 || d1 > n1) stop(sprintf("IV1 인덱스 범위 오류: %02d (유효: 01~%02d)", d1, n1))
  if (d2 < 1 || d2 > n2) stop(sprintf("IV2 인덱스 범위 오류: %02d (유효: 01~%02d)", d2, n2))
  if (d3 < 1 || d3 > n3) stop(sprintf("IV3 인덱스 범위 오류: %02d (유효: 01~%02d)", d3, n3))
  if (d4 < 1 || d4 > n4) stop(sprintf("IV4 인덱스 범위 오류: %02d (유효: 01~%02d)", d4, n4))

  list(
    iv1_score_gap  = IV1_LEVELS[d1],
    iv2_b_interval = IV2_LEVELS[d2],
    iv3_b_mean     = IV3_LEVELS[d3],
    iv4_theta_dist = IV4_LEVELS[d4]
  )
}

# ── 문항 모수 생성 ────────────────────────────────────────────────────────────
#
# 반환값: N_ITEMS 길이의 리스트. 각 원소는 list(sf, b, a)
#
# 저장: output/true_params/YYYYMMDD_cond[C]_rep[R]_seed[S]_item_params.csv
#   열: item, a, b1, b2, b3, b4, sf1, sf2, sf3, sf4, sf5
#
generate_item_params <- function(cond_params, rep_id, seed, cond_code) {

  set.seed(seed)

  d    <- cond_params$iv2_b_interval
  mu   <- cond_params$iv3_b_mean
  sf4  <- cond_params$iv1_score_gap   # 채점함수 네 번째 값

  # ─ 문항 1 (조작 문항) ──────────────────────────────────────────────────────
  # 채점함수: 네 번째 값만 IV1에 따라 변화, 마지막은 4로 고정
  item1_sf <- c(0, 1, 2, sf4, 4)
  item1_b  <- mu + d * c(-1.5, -0.5, 0.5, 1.5)

  # ─ 문항 2-20 (고정 구조) ───────────────────────────────────────────────────
  item_means <- rnorm(N_ITEMS - 1, mean = 0, sd = 1)

  items      <- vector("list", N_ITEMS)
  items[[1]] <- list(sf = item1_sf, b = item1_b, a = DISCRIM)

  for (i in 2:N_ITEMS) {
    m_i        <- item_means[i - 1]
    items[[i]] <- list(
      sf = c(0, 1, 2, 3, 4),
      b  = m_i + 1.5 * c(-1.5, -0.5, 0.5, 1.5),
      a  = DISCRIM
    )
  }

  # ─ CSV 저장 ────────────────────────────────────────────────────────────────
  date_str <- format(Sys.Date(), "%Y%m%d")
  fname <- file.path(
    "output", "true_params",
    sprintf("%s_cond%s_rep%04d_seed%d_item_params.csv",
            date_str, cond_code, rep_id, seed)
  )

  df <- data.frame(
    item = seq_len(N_ITEMS),
    a    = sapply(items, `[[`, "a"),
    b1   = sapply(items, function(x) x$b[1]),
    b2   = sapply(items, function(x) x$b[2]),
    b3   = sapply(items, function(x) x$b[3]),
    b4   = sapply(items, function(x) x$b[4]),
    sf1  = sapply(items, function(x) x$sf[1]),
    sf2  = sapply(items, function(x) x$sf[2]),
    sf3  = sapply(items, function(x) x$sf[3]),
    sf4  = sapply(items, function(x) x$sf[4]),
    sf5  = sapply(items, function(x) x$sf[5])
  )
  write.csv(df, fname, row.names = FALSE)

  items
}
