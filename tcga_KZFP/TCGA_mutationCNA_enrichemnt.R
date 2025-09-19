
# R1: simple check

library(dplyr)
library(purrr)
library(broom)
library(stringr)

# ==== 输入假定 ====
# freq_by_cancer: 每癌种一行，至少包含：
# cohort, group (KZFP_low / KZFP_high), n（该癌种样本数）, 以及若干 *_freq 列（0~1 的频率）

# ---- 统计函数 ----
fisher_from_freq <- function(tbl, col_freq, alt="greater") {
  x <- tbl %>%
    transmute(group,
              pos = round(.data[[col_freq]] * n),
              neg = n - round(.data[[col_freq]] * n)) %>%
    group_by(group) %>%
    summarise(pos=sum(pos), neg=sum(neg), .groups="drop")
  m <- as.matrix(x %>% select(pos, neg)); rownames(m) <- x$group
  m <- m[c("KZFP_low","KZFP_high"), ]
  fisher.test(m, alternative=alt)
}

glm_weighted <- function(tbl, col_freq){
  df <- tbl %>%
    mutate(events = round(.data[[col_freq]] * n),
           non_events = n - events)
  glm(cbind(events, non_events) ~ group,
      data=df, family=binomial(), weights=n)
}

# ---- 识别所有频率列并批量跑 ----
freq_cols <- grep("_freq$", colnames(freq_by_cancer), value=TRUE)

res <- map_dfr(freq_cols, function(colname){
  # Fisher
  f <- fisher_from_freq(freq_by_cancer, colname)
  # Logistic
  m <- glm_weighted(freq_by_cancer, colname)
  tidy_m <- tidy(m) %>% filter(term=="groupKZFP_low")
  # 计算组间“频率差”（low减high，作为直观效应量）
  delta <- freq_by_cancer %>%
    select(group, n, !!sym(colname)) %>%
    summarise(
      low  = weighted.mean(.data[[colname]][group=="KZFP_low"],  w = n[group=="KZFP_low"],  na.rm=TRUE),
      high = weighted.mean(.data[[colname]][group=="KZFP_high"], w = n[group=="KZFP_high"], na.rm=TRUE),
      .groups="drop"
    ) %>%
    mutate(diff = low - high) %>% pull(diff)

  tibble(
    metric   = colname,
    fisher_p = f$p.value,
    logit_OR = exp(tidy_m$estimate),
    logit_p  = tidy_m$p.value,
    diff_low_minus_high = delta
  )
}) %>%
  mutate(
    fisher_FDR = p.adjust(fisher_p, method="BH"),
    logit_FDR  = p.adjust(logit_p,  method="BH")
  )

# 保存
readr::write_tsv(res, "cancer_level_enrichment_results.tsv")
