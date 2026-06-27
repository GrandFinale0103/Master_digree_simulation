# =============================================================================
# 03_response.R
# GPCM 확률 계산 및 가상 응답 생성 함수
# =============================================================================

# ── GPCM 범주 확률 계산 (채점함수 일반화) ────────────────────────────────────
#
# 모형:
#   eta[k] = sum_{j=1}^{k} a * w[j] * (theta - b[j])
#   P(X=k|theta) = exp(eta[k]) / sum_c exp(eta[c])
#
# 여기서 w[j] = sf[j] - sf[j-1] (채점함수 인접 범주 차이)
#
gpcm_sf_prob <- function(theta, a, b, sf) {
  K <- length(sf)
  w <- diff(sf)

  eta    <- numeric(K)
  eta[1] <- 0
  for (k in 2:K) {
    eta[k] <- sum(a * w[seq_len(k - 1)] * (theta - b[seq_len(k - 1)]))
  }

  eta_adj <- eta - max(eta)
  p       <- exp(eta_adj) / sum(exp(eta_adj))
  p
}

# ── 능력모수 표집 ─────────────────────────────────────────────────────────────
sample_theta <- function(n, dist) {
  switch(dist,
    "pos_skew" = sn::rsn(n,
                         xi    = RSN_POS$xi,
                         omega = RSN_POS$omega,
                         alpha = RSN_POS$alpha),
    "normal"   = rnorm(n, mean = 0, sd = 1),
    "neg_skew" = sn::rsn(n,
                         xi    = RSN_NEG$xi,
                         omega = RSN_NEG$omega,
                         alpha = RSN_NEG$alpha),
    "uniform"  = runif(n, min = -4, max = 4),
    stop(sprintf("알 수 없는 분포: '%s'", dist))
  )
}

# ── 빈 범주 확인 ──────────────────────────────────────────────────────────────
#
# 반환: 빈 범주가 있는 문항 목록
#   list(item_index = i, missing_cats = c(...))
#   빈 범주 없으면 빈 리스트 반환
#
find_empty_cats <- function(resp_mat, n_cat) {
  all_cats <- 0:(n_cat - 1)
  empty <- list()
  for (i in seq_len(ncol(resp_mat))) {
    missing <- setdiff(all_cats, unique(resp_mat[, i]))
    if (length(missing) > 0) {
      empty[[length(empty) + 1]] <- list(item = i, missing = missing)
    }
  }
  empty
}

# ── 가상 응답 행렬 생성 ───────────────────────────────────────────────────────
#
# 저장: output/responses/YYYYMMDD_cond[C]_rep[R]_seed[S]_response.csv
#   열: theta (유령 응답자는 NA), item1, ..., item20
#   응답값 범위: 0 ~ N_CAT-1
#
# 빈 범주 처리:
#   빈 범주가 발생하면 유령 응답자를 1명 추가
#   - 빈 범주 문항: 해당 빈 범주를 응답 (복수 빈 범주 시 첫 번째)
#   - 나머지 문항: 0~(N_CAT-1) 무작위 응답
#   빈 범주가 완전히 없어질 때까지 반복 (한 문항에 빈 범주 복수 시 대비)
#
generate_response <- function(item_params, cond_params, rep_id, seed, cond_code) {

  set.seed(seed)

  theta <- sample_theta(N_PERSONS, cond_params$iv4_theta_dist)

  resp <- matrix(NA_integer_, nrow = N_PERSONS, ncol = N_ITEMS)
  colnames(resp) <- paste0("item", seq_len(N_ITEMS))

  for (i in seq_len(N_ITEMS)) {
    ip       <- item_params[[i]]
    prob_mat <- t(vapply(theta, function(th) {
      gpcm_sf_prob(th, ip$a, ip$b, ip$sf)
    }, numeric(N_CAT)))

    resp[, i] <- apply(prob_mat, 1, function(p) {
      sample.int(N_CAT, size = 1, prob = p) - 1L
    })
  }

  # ── 빈 범주 처리: 유령 응답자 추가 ──────────────────────────────────────────
  n_ghost <- 0L

  repeat {
    empty <- find_empty_cats(resp, N_CAT)
    if (length(empty) == 0) break

    # 유령 응답자 응답 벡터 초기화: 모든 문항에 무작위 범주
    ghost_resp <- sample(0:(N_CAT - 1), N_ITEMS, replace = TRUE)

    # 빈 범주가 있는 문항에는 해당 빈 범주 중 첫 번째를 응답
    for (e in empty) {
      ghost_resp[e$item] <- e$missing[1]
    }

    resp  <- rbind(resp, ghost_resp)
    theta <- c(theta, NA_real_)   # 유령 응답자는 theta 없음
    n_ghost <- n_ghost + 1L
  }

  # ── CSV 저장 ─────────────────────────────────────────────────────────────────
  date_str <- format(Sys.Date(), "%Y%m%d")
  fname <- file.path(
    "output", "responses",
    sprintf("%s_cond%s_rep%04d_seed%d_response.csv",
            date_str, cond_code, rep_id, seed)
  )

  df <- cbind(data.frame(theta = theta), as.data.frame(resp))
  write.csv(df, fname, row.names = FALSE)

  list(response = resp, theta = theta, n_ghost = n_ghost)
}
