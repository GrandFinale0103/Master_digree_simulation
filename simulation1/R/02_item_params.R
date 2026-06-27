# =============================================================================
# 02_item_params.R
# 조건 코드 파싱 및 문항 모수 생성 함수
# =============================================================================

# ── 조건 코드 파싱 ────────────────────────────────────────────────────────────
#
# 4자리 문자열 코드를 독립변수 수준값으로 변환
#   예) "1233" → iv1=1, iv2=1.5, iv3=3, iv4="neg_skew"
#
parse_cond_code <- function(code) {
  if (!grepl("^[123][123][123][1234]$", code)) {
    stop(sprintf("유효하지 않은 조건 코드: '%s'\n첫째~셋째 자리는 1-3, 넷째 자리는 1-4", code))
  }
  digits <- as.integer(strsplit(code, "")[[1]])
  list(
    iv1_score_gap  = IV1_LEVELS[digits[1]],
    iv2_b_interval = IV2_LEVELS[digits[2]],
    iv3_b_mean     = IV3_LEVELS[digits[3]],
    iv4_theta_dist = IV4_LEVELS[digits[4]]
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
  sf4  <- cond_params$iv1_score_gap   # 채점함수 네 번째 값 (3, 3.33, 3.66)

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
