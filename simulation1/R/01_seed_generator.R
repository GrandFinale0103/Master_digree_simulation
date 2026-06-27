# =============================================================================
# 01_seed_generator.R
# 전체 시뮬레이션 시드 테이블 사전 생성 — 최초 1회만 실행
#
# 실행 방법: source("R/01_seed_generator.R")
# 주의: output/seed_table.rds 가 이미 존재하면 덮어쓰지 않음
# =============================================================================

source("R/00_setup.R")

SEED_TABLE_PATH <- "output/seed_table.rds"

if (file.exists(SEED_TABLE_PATH)) {
  cat("seed_table.rds 이미 존재함 — 새로 생성하지 않습니다.\n")
  cat("강제로 재생성하려면 파일을 삭제 후 다시 실행하세요.\n")
} else {
  MASTER_SEED <- 20260625
  set.seed(MASTER_SEED)

  N_COND <- length(ALL_COND_CODES)  # 108
  N_REP  <- 1000

  seed_table <- data.frame(
    cond_code = rep(ALL_COND_CODES, each = N_REP),
    rep_id    = rep(seq_len(N_REP), times = N_COND),
    seed      = sample.int(.Machine$integer.max, N_COND * N_REP, replace = FALSE),
    stringsAsFactors = FALSE
  )

  saveRDS(seed_table, SEED_TABLE_PATH)
  cat(sprintf("시드 테이블 생성 완료: %d 행 (%d 조건 × %d 반복)\n",
              nrow(seed_table), N_COND, N_REP))
  cat(sprintf("마스터 시드: %d\n", MASTER_SEED))
  cat(sprintf("저장 경로: %s\n", SEED_TABLE_PATH))
}
