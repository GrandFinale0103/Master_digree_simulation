# =============================================================================
# 12_analysis_nrm_test.R  (simulation2)
# 11_nrm_estimation.R 결과를 모아 반복별 표와 조건별 집계표를 생성
#
# 반복별 결과 (nrm_test_rep.csv):
#   NRM 경계모수 b1~b4, b4−b3, SE, z, p-value(양측), 기각 여부, 전치 여부
#
# 조건별 집계 (nrm_test_cond.csv):
#   n_converged          — 수렴 성공 반복 수
#   n/prop_transposed_34 — b3 ≥ b4 (범주3-4, 4-5 경계 전치) 비율   [분모: 수렴]
#   n/prop_transposed_any— b1~b4 중 하나라도 전치된 비율            [분모: 수렴]
#   n_tested             — 수렴 + SE 계산 성공 반복 수
#   n/prop_reject        — |z| > 1.96 로 H0(b4 = b3) 기각 비율     [분모: n_tested]
#   경계모수 상태(b3 vs b4) × 검정 결과 교차 집계              [모두 n_tested 기준]
#     n_ordered            — 서열화 상태(b3 < b4) 반복 수
#     n_ordered_reject     — 서열화 상태 + H0 기각  (서열이 유의함)
#     n_ordered_not_reject — 서열화 상태 + 기각 못 함
#     prop_reject_ordered  — 서열화 상태 중 기각 비율 (= n_ordered_reject / n_ordered)
#     n_trans              — 전치 상태(b3 ≥ b4) 반복 수
#     n_trans_reject       — 전치 상태 + H0 기각  (전치가 유의함)
#     n_trans_not_reject   — 전치 상태 + 기각 못 함
#     prop_reject_trans    — 전치 상태 중 기각 비율 (= n_trans_reject / n_trans)
#     prop_ordered_reject_of_tested / prop_trans_reject_of_tested — 검정 전체 대비 비율
#   n/prop_sf_rev_34     — 채점함수 전치 (ak3 ≤ ak2 또는 ak4 ≤ ak3: b3·b4 관련) [분모: 수렴]
#   n/prop_sf_transposed_any — ak0~ak4 중 하나라도 전치                  [분모: 수렴]
#
# 출력 폴더: output/nrm/analysis/
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

NRM_DIR <- "output/nrm/estimated_params"
OUT_DIR <- "output/nrm/analysis"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

IV1_MAP <- c("01"=3,"02"=3.2,"03"=3.33,"04"=3.4,"05"=3.6,"06"=3.66,"07"=3.8)
IV2_MAP <- c("01"=1.4,"02"=1.25,"03"=1.2,"04"=1,"05"=0.8,
             "06"=0.6,"07"=0.5,"08"=0.4,"09"=0.25,"10"=0.2)
IV3_MAP <- c("01"=0,"02"=0.5,"03"=1,"04"=1.5)

files <- list.files(NRM_DIR, pattern = "_nrm\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("NRM 결과가 없습니다. 먼저 R/11_nrm_estimation.R 을 실행하세요.")

read_fn <- if (requireNamespace("data.table", quietly = TRUE)) {
  function(p) data.table::fread(p, data.table = FALSE, showProgress = FALSE,
                                colClasses = list(character = c("cond_code", "seed", "error")))
} else {
  function(p) read.csv(p, stringsAsFactors = FALSE,
                       colClasses = c(cond_code = "character", seed = "character", error = "character"))
}
message(sprintf("NRM 결과 %d개 로딩 중...", length(files)))
rep_df <- bind_rows(lapply(files, read_fn))

to_lgl <- function(x) as.logical(toupper(as.character(x)))
rep_df <- rep_df %>%
  mutate(
    across(c(success, converged, se_ok, reject_h0, transposed_34, transposed_any), to_lgl),
    iv1_sf4        = unname(IV1_MAP[substr(cond_code, 1, 2)]),
    iv2_b_interval = unname(IV2_MAP[substr(cond_code, 3, 4)]),
    iv3_b_mean     = unname(IV3_MAP[substr(cond_code, 5, 6)]),
    b3_true = iv3_b_mean + iv2_b_interval * 0.5,
    b4_true = iv3_b_mean + iv2_b_interval * 1.5,
    # 채점함수 전치 (ak 값에서 다시 계산 → 이전 버전 결과 파일도 처리 가능)
    sf_rev_1 = ak1 <= ak0,
    sf_rev_2 = ak2 <= ak1,
    sf_rev_3 = ak3 <= ak2,
    sf_rev_4 = ak4 <= ak3,
    sf_transposed_any = sf_rev_1 | sf_rev_2 | sf_rev_3 | sf_rev_4
  ) %>%
  select(cond_code, rep_id, seed, iv1_sf4, iv2_b_interval, iv3_b_mean,
         success, converged, se_ok,
         a1, ak0, ak1, ak2, ak3, ak4, d0, d1, d2, d3, d4,
         b1, b2, b3, b4, b3_true, b4_true,
         diff_b4_b3, se_diff, z, p_value, reject_h0,
         transposed_34, transposed_any,
         sf_rev_1, sf_rev_2, sf_rev_3, sf_rev_4, sf_transposed_any, error) %>%
  arrange(cond_code, rep_id)

# 양측 검정 판정 (z에서 다시 계산 → 이전 단측 결과 파일도 그대로 사용 가능)
Z_CRIT <- 1.96
rep_df <- rep_df %>%
  mutate(p_value   = 2 * pnorm(abs(z), lower.tail = FALSE),
         reject_h0 = abs(z) > Z_CRIT)

write.csv(rep_df, file.path(OUT_DIR, "nrm_test_rep.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s/nrm_test_rep.csv (%d행)", OUT_DIR, nrow(rep_df)))

summarise_block <- function(df) {
  df %>% summarise(
    n_total             = n(),
    n_success           = sum(success, na.rm = TRUE),
    n_converged         = sum(converged %in% TRUE),
    n_transposed_34     = sum(converged %in% TRUE & transposed_34 %in% TRUE),
    prop_transposed_34  = n_transposed_34 / n_converged,
    n_transposed_any    = sum(converged %in% TRUE & transposed_any %in% TRUE),
    prop_transposed_any = n_transposed_any / n_converged,
    n_sf_rev_34         = sum(converged %in% TRUE & (sf_rev_3 | sf_rev_4) %in% TRUE),
    prop_sf_rev_34      = n_sf_rev_34 / n_converged,
    n_sf_transposed_any = sum(converged %in% TRUE & sf_transposed_any %in% TRUE),
    prop_sf_transposed_any = n_sf_transposed_any / n_converged,
    n_tested            = sum(converged %in% TRUE & se_ok %in% TRUE),
    n_reject            = sum(converged %in% TRUE & reject_h0 %in% TRUE),
    prop_reject         = n_reject / n_tested,
    n_ordered            = sum(converged %in% TRUE & se_ok %in% TRUE & transposed_34 %in% FALSE),
    n_ordered_reject     = sum(converged %in% TRUE & se_ok %in% TRUE & transposed_34 %in% FALSE & reject_h0 %in% TRUE),
    n_ordered_not_reject = n_ordered - n_ordered_reject,
    prop_reject_ordered  = n_ordered_reject / n_ordered,
    n_trans              = sum(converged %in% TRUE & se_ok %in% TRUE & transposed_34 %in% TRUE),
    n_trans_reject       = sum(converged %in% TRUE & se_ok %in% TRUE & transposed_34 %in% TRUE & reject_h0 %in% TRUE),
    n_trans_not_reject   = n_trans - n_trans_reject,
    prop_reject_trans    = n_trans_reject / n_trans,
    prop_ordered_reject_of_tested = n_ordered_reject / n_tested,
    prop_trans_reject_of_tested   = n_trans_reject / n_tested,
    mean_z              = mean(z[converged %in% TRUE], na.rm = TRUE),
    mean_b3             = mean(b3[converged %in% TRUE], na.rm = TRUE),
    mean_b4             = mean(b4[converged %in% TRUE], na.rm = TRUE),
    .groups = "drop"
  )
}

cond_df <- rep_df %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean) %>%
  summarise_block() %>%
  mutate(b3_true = iv3_b_mean + iv2_b_interval * 0.5,
         b4_true = iv3_b_mean + iv2_b_interval * 1.5,
         .after = iv3_b_mean)

write.csv(cond_df, file.path(OUT_DIR, "nrm_test_cond.csv"), row.names = FALSE)
message(sprintf("저장 완료: %s/nrm_test_cond.csv (%d 조건)", OUT_DIR, nrow(cond_df)))

overall <- summarise_block(rep_df)
sink(file.path(OUT_DIR, "nrm_test_overall.txt"))
tryCatch({
  cat("NRM 경계모수 b3 vs b4 양측 검정\n")
  cat("H0: b4 = b3   HA: b4 != b3   기각: |z| = |(b4 - b3)/SE| > 1.96\n\n")
  print(as.data.frame(t(overall)), digits = 4)
}, finally = sink())
message(sprintf("저장 완료: %s/nrm_test_overall.txt", OUT_DIR))
