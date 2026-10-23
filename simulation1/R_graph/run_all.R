# =============================================================================
# run_all.R  (simulation1/R_graph)
# 그래프·표를 순서대로 모두 만듦
#
# 사용법: 작업 디렉토리를 simulation1/ 로 맞춘 뒤
#           source("R_graph/run_all.R", encoding = "UTF-8")
#
# COLLECT = FALSE 로 두면 01_collect_data.R (결과 파일 읽기, 가장 오래 걸림)를 건너뛰고
# 이미 저장된 <GRAPH_DIR>/data/*.rds 를 사용 (GRAPH_DIR 은 00_config.R 에서 정함)
# =============================================================================

COLLECT <- TRUE

scripts <- c(
  if (COLLECT) "01_collect_data.R",
  "02_fig_concept.R",       # 그림 1~5   연구 목적
  "03_fig_design.R",        # 그림 6     방법
  "04_fig_transposition.R", # 그림 7·8·12 전치 비율
  "05_fig_error.R",         # 그림 9·13  편향·RMSE
  "06_fig_mechanism.R",     # 그림 10·11 편향 원리
  "07_fig_nrm.R",           # 그림 14·15 NRM 전치 원인
  "08_fig_summary.R"        # 그림 16    요약 히트맵
)

failed <- character(0)
for (s in scripts) {
  message("\n==== ", s, " ====")
  # 스크립트마다 새 환경에서 실행 (변수 충돌 방지)
  tryCatch(source(file.path("R_graph", s), local = new.env(parent = globalenv()),
                  encoding = "UTF-8"),
           error = function(e) {
             message("  [오류] ", s, ": ", conditionMessage(e))
             failed <<- c(failed, sprintf("%s: %s", s, conditionMessage(e)))
           })
}
if (length(failed)) {
  message("\n[주의] 실패한 스크립트가 있습니다. 이 스크립트의 그림은 새로 그려지지 않았고,")
  message("       output/figures 에 예전 그림이 남아 있을 수 있습니다:")
  for (f in failed) message("  - ", f)
} else {
  .cfg <- new.env(); sys.source("R_graph/00_config.R", envir = .cfg, keep.source = FALSE)
  message(sprintf("\n완료: 그림 %s/, 표 %s/", .cfg$FIG_DIR, .cfg$TAB_DIR))
}
