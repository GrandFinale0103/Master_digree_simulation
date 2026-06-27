# =============================================================================
# 04_estimation.R
# mirt 패키지를 이용한 1차원 GPCM 모수 추정 함수
# =============================================================================

# ── GPCM 모수 추정 ────────────────────────────────────────────────────────────
#
# 추정 모형: 1차원 GPCM (표준 채점함수 [0,1,2,3,4] 가정)
#   → 실제 데이터가 비등간격 채점함수로 생성되었더라도
#     등간격을 가정하고 추정 (연구 목적: 이 불일치의 영향 탐색)
#
# 저장: output/estimated_params/YYYYMMDD_cond[C]_rep[R]_seed[S]_est_params.csv
#   성공 시: 20행 × (item, converged, a, b1, b2, b3, b4)
#   실패 시: 1행  × (success=FALSE, error=에러메시지)
#
estimate_params <- function(response_matrix, rep_id, seed, cond_code) {

  date_str <- format(Sys.Date(), "%Y%m%d")
  fname <- file.path(
    "output", "estimated_params",
    sprintf("%s_cond%s_rep%04d_seed%d_est_params.csv",
            date_str, cond_code, rep_id, seed)
  )

  result <- tryCatch({
    resp_df <- as.data.frame(response_matrix + 1L)

    mod <- mirt::mirt(
      data      = resp_df,
      model     = 1,
      itemtype  = "gpcm",
      verbose   = FALSE,
      technical = list(NCYCLES = 3000)
    )

    converged <- mirt::extract.mirt(mod, "converged")
    coef_irt  <- mirt::coef(mod, simplify = FALSE, IRTpars = TRUE)

    # GroupPars 제외하고 문항 계수만 추출
    item_names <- names(coef_irt)[names(coef_irt) != "GroupPars"]
    df <- do.call(rbind, lapply(item_names, function(nm) {
      co <- coef_irt[[nm]]
      data.frame(
        item      = nm,
        converged = converged,
        a         = co[1, "a"],
        b1        = co[1, "b1"],
        b2        = co[1, "b2"],
        b3        = co[1, "b3"],
        b4        = co[1, "b4"],
        stringsAsFactors = FALSE
      )
    }))
    write.csv(df, fname, row.names = FALSE)

    list(success = TRUE, coef = coef_irt, convergence = converged)

  }, error = function(e) {
    df <- data.frame(success = FALSE, error = conditionMessage(e),
                     stringsAsFactors = FALSE)
    write.csv(df, fname, row.names = FALSE)
    list(success = FALSE, error = conditionMessage(e))
  })

  result
}
