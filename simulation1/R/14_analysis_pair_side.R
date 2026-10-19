# =============================================================================
# 14_analysis_pair_side.R  (simulation1, 정규분포 조건만)
# 문항 1에서 인접한 두 범주 중 하나를 응답한 사람들이
# 참 경계모수(true b) 기준으로 아래쪽(θ < b)과 위쪽(θ ≥ b) 중 어디에 더 많은지 집계
#
# ── 범주와 경계모수 ───────────────────────────────────────────────────────────
#   응답 파일의 응답값은 0~4. 범주를 1~5로 세면 응답값 + 1.
#   경계 b_k 는 범주 k 와 k+1 (응답값 k−1 과 k) 사이의 경계
#     예) 범주 4·5 (응답값 3·4) ↔ b4
#   참 경계모수: b = IV3 + IV2 × (−1.5, −0.5, 0.5, 1.5)  (문항 1은 반복마다 같음)
#
# ── 집계 (반복마다, 유령 응답자(θ = NA) 제외) ─────────────────────────────────
#   n_pair      : 범주 k 또는 k+1 을 응답한 사람 수
#   n_below     : 그중 θ <  b_k 인 사람 수
#   n_above     : 그중 θ >= b_k 인 사람 수
#   prop_below  : n_below / n_pair
#   mean_theta  : 그 사람들의 평균 θ
#   조건별: 반복을 합친 비율(prop_below_pooled = Σn_below / Σn_pair),
#           반복별 비율의 평균·표준편차, 평균 θ
#
# ── 이론값 (표본이 무한할 때) ─────────────────────────────────────────────────
#   θ ~ N(0, 1) 과 생성 모형(채점함수 0, 1, 2, s4, 4, a = 1)으로 적분해 계산
#   prop_below_theory = P(θ < b_k | 범주 k 또는 k+1 응답)
#   mean_theta_theory = E(θ | 범주 k 또는 k+1 응답)
#   → 시뮬레이션 집계가 이 값 근처에 있어야 함 (집계 검증용)
#
# 입력: output/responses/*_cond???2_rep*_response.csv  (조건 코드 넷째 자리 2 = 정규분포)
# 출력: output/analysis_normal/pair_side/
#   pair_side_rep.csv          — 반복 × 경계별 원자료
#   pair_side_cond.csv         — 조건 × 경계별 요약 (이론값 포함)
#   pair_side_summary.txt      — 경계별 요약 표 (콘솔 출력과 같음)
#   pair_side_plot_b{k}.png    — 경계별 아래쪽 비율 (행 IV1, 열 IV2, x축 IV3)
#   pair_side_bar_b{k}.png     — 경계별 아래쪽/위쪽 구성 누적 막대 (같은 배치)
# =============================================================================

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                          사용자 설정 영역                                 ║
# ╠══════════════════════════════════════════════════════════════════════════╣
PAIRS   <- 1:4      # 집계할 경계 (4 = 범주 4·5). b4만 보려면 PAIRS <- 4
ITEM    <- "item1"  # 집계할 문항 (조작 문항)
N_CORES <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
# ╚══════════════════════════════════════════════════════════════════════════╝

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(parallel)
})

RESP_DIR <- "output/responses"
OUT_DIR  <- "output/analysis_normal/pair_side"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(OUT_DIR, f)

IV1_MAP  <- c("1" = 3.00, "2" = 3.33, "3" = 3.66)
IV2_MAP  <- c("1" = 1.5,  "2" = 1.0,  "3" = 0.5)
IV3_MAP  <- c("1" = 0.0,  "2" = 0.5,  "3" = 1.0)
B_OFFSET <- c(-1.5, -0.5, 0.5, 1.5)

true_b <- function(cond_code) {
  unname(IV3_MAP[substr(cond_code, 3, 3)]) +
    unname(IV2_MAP[substr(cond_code, 2, 2)]) * B_OFFSET
}

# =============================================================================
# 1. 응답 파일 목록 (정규분포 조건만)
# =============================================================================
files <- list.files(RESP_DIR, pattern = "_cond[0-9]{3}2_rep[0-9]+_.*_response\\.csv$",
                    full.names = TRUE)
if (length(files) == 0) stop(sprintf("%s 에 정규분포 조건 응답 파일이 없습니다.", RESP_DIR))
message(sprintf("응답 파일 %d 개 (정규분포 조건 %d 개)", length(files),
                length(unique(sub(".*_cond([0-9]+)_.*", "\\1", basename(files))))))

# =============================================================================
# 2. 반복별 집계
# =============================================================================
conds  <- unique(sub(".*_cond([0-9]+)_.*", "\\1", basename(files)))
b_list <- setNames(lapply(conds, true_b), conds)   # 조건별 참 경계모수 (병렬 작업자에 그대로 전달)

count_one <- function(path, pairs, item, b_list) {
  b_name <- basename(path)
  cond   <- sub(".*_cond([0-9]+)_.*", "\\1", b_name)
  rep_id <- as.integer(sub(".*_rep([0-9]+)_.*", "\\1", b_name))
  d <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, select = c("theta", item), data.table = FALSE)
  } else {
    read.csv(path)[, c("theta", item)]
  }
  d <- d[!is.na(d$theta), ]                     # 유령 응답자 제외
  b <- b_list[[cond]]
  do.call(rbind, lapply(pairs, function(k) {
    th <- d$theta[d[[item]] %in% c(k - 1, k)]   # 응답값 k−1, k = 범주 k, k+1
    data.frame(cond_code = cond, rep_id = rep_id, boundary = k, b_true = b[k],
               n_total = nrow(d), n_pair = length(th),
               n_below = sum(th < b[k]), n_above = sum(th >= b[k]),
               mean_theta = if (length(th)) mean(th) else NA_real_)
  }))
}

message(sprintf("집계 중 (%d 코어)...", N_CORES))
if (N_CORES > 1) {
  cl <- makeCluster(N_CORES)
  rep_list <- parLapply(cl, files, count_one, pairs = PAIRS, item = ITEM, b_list = b_list)
  stopCluster(cl)
} else {
  rep_list <- lapply(files, count_one, pairs = PAIRS, item = ITEM, b_list = b_list)
}
rep_df <- bind_rows(rep_list) %>%
  mutate(prop_below = ifelse(n_pair > 0, n_below / n_pair, NA_real_),
         iv1_sf4        = unname(IV1_MAP[substr(cond_code, 1, 1)]),
         iv2_b_interval = unname(IV2_MAP[substr(cond_code, 2, 2)]),
         iv3_b_mean     = unname(IV3_MAP[substr(cond_code, 3, 3)])) %>%
  arrange(cond_code, boundary, rep_id)
write.csv(rep_df, out("pair_side_rep.csv"), row.names = FALSE)

# =============================================================================
# 3. 이론값 (θ ~ N(0, 1), 생성 모형으로 적분)
# =============================================================================
QUAD   <- seq(-6, 6, length.out = 1201)
QUAD_W <- dnorm(QUAD) / sum(dnorm(QUAD))

true_probs <- function(theta, sf, b) {          # 03_response.R 의 gpcm_sf_prob 와 같은 식 (a = 1)
  w   <- diff(sf)
  eta <- outer(theta, c(0, cumsum(w))) - matrix(c(0, cumsum(w * b)), length(theta), 5, byrow = TRUE)
  eta <- eta - apply(eta, 1, max)
  p   <- exp(eta); p / rowSums(p)
}

theory <- expand.grid(cond_code = unique(rep_df$cond_code), boundary = PAIRS,
                      stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(th = list({
    b  <- true_b(cond_code)
    sf <- c(0, 1, 2, unname(IV1_MAP[substr(cond_code, 1, 1)]), 4)
    P  <- true_probs(QUAD, sf, b)
    wt <- QUAD_W * (P[, boundary] + P[, boundary + 1])   # 범주 k, k+1 응답 확률
    c(prop = sum(wt[QUAD < b[boundary]]) / sum(wt), mean = sum(wt * QUAD) / sum(wt))
  }),
  prop_below_theory = th[["prop"]], mean_theta_theory = th[["mean"]]) %>%
  ungroup() %>% select(-th)

# =============================================================================
# 4. 조건별 요약
# =============================================================================
cond_df <- rep_df %>%
  group_by(cond_code, iv1_sf4, iv2_b_interval, iv3_b_mean, boundary, b_true) %>%
  summarise(n_rep             = n(),
            mean_n_pair       = mean(n_pair),
            mean_n_below      = mean(n_below),
            mean_n_above      = mean(n_above),
            prop_below_pooled = sum(n_below) / sum(n_pair),
            prop_below_mean   = mean(prop_below, na.rm = TRUE),
            prop_below_sd     = sd(prop_below, na.rm = TRUE),
            mean_theta        = weighted.mean(mean_theta, n_pair, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(prop_above_pooled = 1 - prop_below_pooled,
         more_on           = ifelse(prop_below_pooled > 0.5, "below", "above")) %>%
  left_join(theory, by = c("cond_code", "boundary")) %>%
  arrange(boundary, iv1_sf4, desc(iv2_b_interval), iv3_b_mean)
write.csv(cond_df, out("pair_side_cond.csv"), row.names = FALSE)

# =============================================================================
# 5. 요약 출력
# =============================================================================
sink_file <- out("pair_side_summary.txt")
sink(sink_file, split = TRUE)
cat("범주 k·k+1 응답자 중 참 경계 b_k 보다 아래(θ < b_k)에 있는 비율 — 정규분포 조건\n")
cat("prop_below: 반복 합산 비율 | theory: θ ~ N(0,1) 이론값 | mean_θ: 해당 응답자 평균 θ\n")
for (k in PAIRS) {
  cat(sprintf("\n== b%d (범주 %d·%d) ==\n", k, k, k + 1))
  tab <- cond_df %>% filter(boundary == k) %>%
    transmute(cond = cond_code, s4 = iv1_sf4, interval = iv2_b_interval, mean = iv3_b_mean,
              b_true = round(b_true, 2), n_pair = round(mean_n_pair, 1),
              prop_below = round(prop_below_pooled, 3), theory = round(prop_below_theory, 3),
              mean_theta = round(mean_theta, 2), more_on)
  print(as.data.frame(tab), row.names = FALSE)
  cat(sprintf("  아래쪽이 더 많은 조건: %d / %d\n", sum(tab$more_on == "below"), nrow(tab)))
}
sink()

# =============================================================================
# 6. 그래프 (경계별 파일, 행 IV1 · 열 IV2 · x축 IV3)
# =============================================================================
plot_df <- cond_df %>%
  mutate(iv1_lab = factor(sprintf("s4 = %.2f", iv1_sf4)),
         iv2_lab = factor(sprintf("interval = %.1f", iv2_b_interval),
                          levels = sprintf("interval = %.1f", sort(unique(iv2_b_interval), decreasing = TRUE))))

for (k in PAIRS) {
  d <- plot_df %>% filter(boundary == k)

  p <- ggplot(d, aes(iv3_b_mean)) +
    geom_hline(yintercept = 0.5, linetype = "dotted", colour = "grey40") +
    geom_errorbar(aes(ymin = pmax(prop_below_mean - prop_below_sd, 0),
                      ymax = pmin(prop_below_mean + prop_below_sd, 1)),
                  width = 0.06, colour = "grey55") +
    geom_line(aes(y = prop_below_pooled, colour = "simulation")) +
    geom_point(aes(y = prop_below_pooled, colour = "simulation"), size = 2.4) +
    geom_point(aes(y = prop_below_theory, colour = "theory (N -> inf)"), shape = 4, size = 3, stroke = 1.1) +
    scale_colour_manual(NULL, values = c("simulation" = "#2a6fb0", "theory (N -> inf)" = "#c0392b")) +
    scale_x_continuous(breaks = unname(IV3_MAP)) +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(iv1_lab ~ iv2_lab) +
    labs(title = sprintf("Item 1: respondents in categories %d or %d -- share below true b%d", k, k + 1, k),
         subtitle = "point = pooled share over reps, bar = +-1 SD across reps, dotted = 0.5",
         x = "IV3 (mean of boundaries)", y = sprintf("P(theta < b%d | category %d or %d)", k, k, k + 1)) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
  ggsave(out(sprintf("pair_side_plot_b%d.png", k)), p, width = 10, height = 8, dpi = 150)

  bar_df <- bind_rows(
    d %>% transmute(iv1_lab, iv2_lab, iv3_b_mean, side = "below b (theta < b)", n = mean_n_below),
    d %>% transmute(iv1_lab, iv2_lab, iv3_b_mean, side = "above b (theta >= b)", n = mean_n_above)) %>%
    mutate(side = factor(side, levels = c("above b (theta >= b)", "below b (theta < b)")))
  pb <- ggplot(bar_df, aes(factor(iv3_b_mean), n, fill = side)) +
    geom_col(width = 0.7) +
    scale_fill_manual(NULL, values = c("below b (theta < b)" = "#2a6fb0", "above b (theta >= b)" = "#e3a33b")) +
    facet_grid(iv1_lab ~ iv2_lab) +
    labs(title = sprintf("Item 1: mean number of respondents in categories %d or %d, split at true b%d", k, k + 1, k),
         x = "IV3 (mean of boundaries)", y = "respondents per replication") +
    theme_bw(base_size = 13) + theme(legend.position = "bottom")
  ggsave(out(sprintf("pair_side_bar_b%d.png", k)), pb, width = 10, height = 8, dpi = 150)
}

message(sprintf("완료: %s", OUT_DIR))
