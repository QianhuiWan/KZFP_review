# ---- 安装依赖（如未装）----
# install.packages(c("tidyverse","data.table"))
# BiocManager::install(c("TCGAbiolinks","maftools"))
# install.packages("msigdbr")

library(tidyverse)
library(data.table)
library(TCGAbiolinks)
library(msigdbr)

# ========= A. 读入你的癌种分组 =========
# 你的文件：两列 cohort, group（如 TCGA-BRCA, KZFP_low / KZFP_high；名称随意但分两组）
group_df <- readr::read_tsv("cancer_groups.tsv") %>%
  mutate(group = as.character(group))

# 可选：统一组名
low_name  <- unique(group_df$group)[1]   # 根据你的表
high_name <- unique(group_df$group)[2]

# ========= B. 基因集 =========
ifn_genes_file <- "ifn_genes.txt"
if (file.exists(ifn_genes_file)) {
  ifn_genes <- toupper(readr::read_lines(ifn_genes_file))
} else {
  ifn_genes <- toupper(c(
    "IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1","IL10RB",
    "JAK1","JAK2","TYK2","STAT1","STAT2","IRF9","IRF1","IRF7","IRF8"
  ))
}
nonsyn <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
            "Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
            "Splice_Site","Translation_Start_Site")
norm_barcode <- function(x) substr(x, 1, 15)

# ========= C. 突变频率（TCGAbiolinks MAF）=========
calc_mut_freq_one <- function(project){
  message("MAF: ", project)
  tumor <- sub("^TCGA-", "", project)
  maf_df <- tryCatch(GDCquery_Maf(tumor = tumor, pipelines = "mutect2"),
                     error=function(e) NULL)
  if (is.null(maf_df) || nrow(maf_df)==0) return(NULL)
  maf_ns <- maf_df %>%
    filter(Variant_Classification %in% nonsyn) %>%
    mutate(Hugo_Symbol = toupper(Hugo_Symbol),
           barcode = norm_barcode(Tumor_Sample_Barcode))
  if (nrow(maf_ns)==0) return(NULL)

  # TP53
  tp53_bar <- maf_ns %>% filter(Hugo_Symbol=="TP53") %>% distinct(barcode)
  # IFN 集合任一基因
  ifn_bar  <- maf_ns %>% filter(Hugo_Symbol %in% ifn_genes) %>% distinct(barcode)

  # denominator：该项目所有突变样本条码集合（也可换成临床样本全集合）
  all_bar  <- maf_ns %>% distinct(barcode) %>% pull()

  tibble(
    cohort       = project,
    n            = length(all_bar),
    tp53_mut_n   = length(intersect(all_bar, tp53_bar$barcode)),
    ifn_mut_n    = length(intersect(all_bar, ifn_bar$barcode))
  ) %>%
    mutate(tp53_mut_freq = tp53_mut_n / n,
           ifn_mut_freq  = ifn_mut_n  / n)
}

mut_list <- lapply(unique(group_df$cohort), calc_mut_freq_one)
mut_freq <- bind_rows(mut_list)

# ========= D. CNV 频率（两种方案二选一）=========

## 方案 A：用 UCSC Xena 的 GISTIC 阈值矩阵（推荐省事）
# （如果允许混用数据源的话）
# library(UCSCXenaTools)
# xcnv <- XenaGenerate() %>% XenaFilter(filterDatasets = "TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes")
# cnv_file <- XenaDownload(xcnv) %>% XenaPrepare()
# cnv_mat <- as.data.frame(cnv_file) %>% column_to_rownames("gene") %>% as.matrix()
# colnames(cnv_mat) <- norm_barcode(colnames(cnv_mat))

# get_cnv_freq <- function(project, gene_set, loss_thr=-1, amp_thr=1){
#   bar_in_project <- colnames(cnv_mat)[grepl(project, colnames(cnv_mat))]  # 如果列名含项目标识，否则需映射
#   if (length(bar_in_project)==0) return(NULL)
#   gs <- intersect(gene_set, rownames(cnv_mat))
#   if (length(gs)==0) return(NULL)
#   sub <- cnv_mat[gs, bar_in_project, drop=FALSE]
#   loss_any <- colSums(sub <= loss_thr, na.rm=TRUE) > 0
#   amp_any  <- colSums(sub >= amp_thr, na.rm=TRUE) > 0
#   tibble(cohort=project,
#          n=length(bar_in_project),
#          loss_n = sum(loss_any), loss_freq = mean(loss_any),
#          amp_n  = sum(amp_any),  amp_freq  = mean(amp_any))
# }

## 方案 B：全用 TCGAbiolinks，从 CNV 段阈值化（纯 GDC 路线）
# 需要你自备 gene 坐标（GRCh38），这里给出一个最简方案：直接统计 seg.mean 阈值命中率（不映射到基因）
# 若要“按基因集合”（IFN/TP53）严格统计，需先把段落到基因，可在下一步我再给你补全坐标版。

# ========= D. CNV 频率（精准到 TP53/IFN 基因集合；可保留宽口径作对照）=========
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(biomaRt)
})

# --- D0) 基因坐标（GRCh38）构建工具 ---
# 返回 GRanges，mcols 包含 symbol
build_gene_coords <- function(genes, genome_version = "GRCh38"){
  genes <- unique(toupper(genes))
  mart <- biomaRt::useEnsembl(biomart = "genes",
                              dataset = "hsapiens_gene_ensembl",
                              version = NULL)  # 当前 GRCh38
  tbl <- biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name",
                                       "start_position","end_position"),
                        filters = "hgnc_symbol",
                        values  = genes,
                        mart    = mart)
  # 只保留常染色体与 X/Y/MT
  tbl <- tbl %>% dplyr::filter(chromosome_name %in% c(as.character(1:22),"X","Y","MT"))
  gr <- GenomicRanges::GRanges(
    seqnames = tbl$chromosome_name,
    ranges   = IRanges::IRanges(start = tbl$start_position, end = tbl$end_position),
    symbol   = toupper(tbl$hgnc_symbol)
  )
  gr
}

# --- D1) 段→基因 overlap + 频率统计（核心函数） ---
# cnv_seg: GDCprepare 的 Copy Number Segment data.frame（含 Chromosome/Start/End/Segment_Mean/Sample）
# gene_coords: GRanges（含 mcols$symbol）
# gene_set: 字符向量（基因符号）
# 返回：cohort, n, loss_n, amp_n, loss_freq, amp_freq
get_cnv_freq_for_genes <- function(project, cnv_seg, gene_coords, gene_set,
                                   loss_thr = -0.3, amp_thr = 0.3){
  message("CNV segment → genes [", project, "]  (|set|=", length(gene_set), ")")
  if (is.null(cnv_seg) || nrow(cnv_seg)==0) return(NULL)

  cnv_seg <- as.data.frame(cnv_seg)
  cnv_seg$barcode <- substr(cnv_seg$Sample, 1, 15)

  seg_gr <- GenomicRanges::GRanges(
    seqnames = cnv_seg$Chromosome,
    ranges   = IRanges::IRanges(start = cnv_seg$Start, end = cnv_seg$End),
    seg.mean = cnv_seg$Segment_Mean,
    sample   = cnv_seg$barcode
  )

  gs <- gene_coords[ gene_coords$symbol %in% toupper(gene_set) ]
  if (length(gs) == 0) {
    warning("No gene coords matched for set in ", project)
    return(NULL)
  }

  hits <- GenomicRanges::findOverlaps(gs, seg_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(tibble::tibble(cohort=project, n=0, loss_n=0, amp_n=0,
                          loss_freq=NA_real_, amp_freq=NA_real_))
  }

  df <- dplyr::tibble(
    gene   = gs$symbol[S4Vectors::queryHits(hits)],
    sample = S4Vectors::mcols(seg_gr)$sample[S4Vectors::subjectHits(hits)],
    seg    = S4Vectors::mcols(seg_gr)$seg.mean[S4Vectors::subjectHits(hits)]
  )

  by_sample <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      loss_any = as.integer(any(seg <= loss_thr, na.rm=TRUE)),
      amp_any  = as.integer(any(seg >= amp_thr,  na.rm=TRUE)),
      .groups="drop"
    )

  tibble::tibble(
    cohort    = project,
    n         = nrow(by_sample),
    loss_n    = sum(by_sample$loss_any),
    amp_n     = sum(by_sample$amp_any)
  ) %>%
    dplyr::mutate(loss_freq = ifelse(n>0, loss_n/n, NA_real_),
                  amp_freq  = ifelse(n>0, amp_n/n,  NA_real_))
}

# --- D2) 仍保留原“宽口径”任意段命中版（可做敏感性分析/对照） ---
get_cnv_seg_freq_any <- function(project, loss_thr=-0.3, amp_thr=0.3){
  message("CNV segment (ANY) [", project, "]")
  q.cnv <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Copy Number Variation",
                                  data.type = "Copy Number Segment")
  cnv_seg <- tryCatch({
    TCGAbiolinks::GDCdownload(q.cnv, method="api", files.per.chunk=50)
    TCGAbiolinks::GDCprepare(q.cnv)
  }, error=function(e) NULL)
  if (is.null(cnv_seg) || nrow(cnv_seg)==0) return(NULL)

  cnv_seg <- as.data.frame(cnv_seg)
  cnv_seg$barcode <- substr(cnv_seg$Sample, 1, 15)

  by_sample <- cnv_seg %>%
    dplyr::group_by(barcode) %>%
    dplyr::summarise(
      loss_any = as.integer(any(Segment_Mean <= loss_thr, na.rm=TRUE)),
      amp_any  = as.integer(any(Segment_Mean >= amp_thr,  na.rm=TRUE)),
      .groups="drop"
    )

  tibble::tibble(
    cohort    = project,
    n         = nrow(by_sample),
    loss_n    = sum(by_sample$loss_any),
    amp_n     = sum(by_sample$amp_any)
  ) %>%
    dplyr::mutate(loss_freq = ifelse(n>0, loss_n/n, NA_real_),
                  amp_freq  = ifelse(n>0, amp_n/n,  NA_real_)) %>%
    dplyr::rename(any_loss_freq = loss_freq, any_amp_freq = amp_freq)
}

get_cnv_seg_freq_any <- function(project, loss_thr=-0.3, amp_thr=0.3){
  message("CNV segment: ", project)
  q.cnv <- GDCquery(project = project,
                    data.category = "Copy Number Variation",
                    data.type = "Copy Number Segment")
  cnv_seg <- tryCatch({
    GDCdownload(q.cnv, method="api", files.per.chunk=50)
    GDCprepare(q.cnv)
  }, error=function(e) NULL)
  if (is.null(cnv_seg) || nrow(cnv_seg)==0) return(NULL)

  cnv_seg <- as.data.frame(cnv_seg)
  cnv_seg$barcode <- norm_barcode(cnv_seg$Sample)

  # 任一段达到缺失/扩增阈值即记为该样本“有缺失/有扩增”
  by_sample <- cnv_seg %>%
    group_by(barcode) %>%
    summarise(
      loss_any = as.integer(any(Segment_Mean <= loss_thr, na.rm=TRUE)),
      amp_any  = as.integer(any(Segment_Mean >= amp_thr,  na.rm=TRUE)),
      .groups="drop"
    )

  tibble(
    cohort    = project,
    n         = nrow(by_sample),
    loss_n    = sum(by_sample$loss_any),
    amp_n     = sum(by_sample$amp_any)
  ) %>%
    mutate(loss_freq = loss_n / n,
           amp_freq  = amp_n  / n)
}

# 注意：上面的“按段任何位置即命中”是“宽松 CNV 频率”，不区分基因集合。
# 如果你需要**只看 TP53 / IFN 基因集合的 CNV**，我可以给你“段→基因（GRCh38坐标）”的版本:


# --- D3) 主循环：下载项目 CNV 段 + 计算 TP53 与 IFN 集合的 CNV 频率 ---
# 1) 预先构建需要的基因坐标（只需一次）
tp53_coords <- build_gene_coords("TP53")
ifn_coords  <- build_gene_coords(ifn_genes)

# 2) 按你的分组表里的癌种循环
get_project_cnv_segments <- function(project){
  q.cnv <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Copy Number Variation",
                                  data.type = "Copy Number Segment")
  tryCatch({
    TCGAbiolinks::GDCdownload(q.cnv, method="api", files.per.chunk=50)
    TCGAbiolinks::GDCprepare(q.cnv)
  }, error=function(e) NULL)
}

cnv_precise_list <- lapply(unique(group_df$cohort), function(prj){
  cnv_seg <- get_project_cnv_segments(prj)
  if (is.null(cnv_seg) || nrow(cnv_seg)==0) {
    return(dplyr::bind_rows(
      tibble::tibble(cohort=prj, n=0, loss_n=0, amp_n=0, loss_freq=NA_real_, amp_freq=NA_real_) %>%
        dplyr::rename(tp53_loss_freq=loss_freq, tp53_amp_freq=amp_freq),
      tibble::tibble(cohort=prj, n=0, loss_n=0, amp_n=0, loss_freq=NA_real_, amp_freq=NA_real_) %>%
        dplyr::rename(ifn_loss_freq=loss_freq,  ifn_amp_freq=amp_freq)
    ))
  }
  tp53 <- get_cnv_freq_for_genes(prj, cnv_seg, tp53_coords, gene_set = "TP53")
  ifn  <- get_cnv_freq_for_genes(prj, cnv_seg, ifn_coords,  gene_set = ifn_genes)

  tp53 <- tp53 %>% dplyr::rename(tp53_loss_freq = loss_freq, tp53_amp_freq = amp_freq)
  ifn  <- ifn  %>% dplyr::rename(ifn_loss_freq  = loss_freq, ifn_amp_freq  = amp_freq)

  # 同一个 project 的 n 可能不同（因为集合不同），我们保留各自 n_tp53 / n_ifn，方便后面检查
  tp53 <- tp53 %>% dplyr::rename(n_tp53 = n)
  ifn  <- ifn  %>% dplyr::rename(n_ifn  = n)

  dplyr::full_join(tp53, ifn, by = "cohort")
})

cnv_precise <- dplyr::bind_rows(cnv_precise_list)

# （可选）同时计算“宽口径 ANY”作为敏感性分析对照
cnv_any_list <- lapply(unique(group_df$cohort), get_cnv_seg_freq_any)
cnv_any <- dplyr::bind_rows(cnv_any_list)



cnv_list <- lapply(unique(group_df$cohort), get_cnv_seg_freq_any)
cnv_freq <- bind_rows(cnv_list) %>% 
  rename(any_loss_freq = loss_freq, any_amp_freq = amp_freq)

# ========= E. 合并、按你给的分组做富集 =========
freq_by_cancer <- group_df %>%
  left_join(mut_freq,   by="cohort") %>%              # 突变频率（已算好）
  left_join(cnv_precise, by="cohort") %>%             # 精准 TP53/IFN 基因集合的 CNV 频率
  left_join(cnv_any,    by="cohort")                  # （可选）宽口径 ANY 频率，对照用

readr::write_tsv(freq_by_cancer, "freq_by_cancer.tsv")


readr::write_tsv(freq_by_cancer, "freq_by_cancer.tsv")

fisher_one_sided <- function(tbl, col_freq, alt = "greater"){
  x <- tbl %>% 
    transmute(group, pos = round(.data[[col_freq]] * n), neg = n - round(.data[[col_freq]] * n)) %>%
    group_by(group) %>% summarise(pos=sum(pos), neg=sum(neg), .groups="drop")
  if (nrow(x) < 2) return(NA_real_)
  m <- as.matrix(x[,c("pos","neg")])
  rownames(m) <- x$group
  if (!all(c(low_name, high_name) %in% rownames(m))) return(NA_real_)
  fisher.test(m[c(low_name, high_name], ), alternative = alt)$p.value
}

targets_mut <- c("tp53_mut_freq","ifn_mut_freq")
# 记得把富集检验的指标换成精准列名：
targets_cnv <- c("tp53_loss_freq","tp53_amp_freq","ifn_loss_freq","ifn_amp_freq")

res_mut <- tibble(metric = targets_mut) %>%
  rowwise() %>% mutate(p = fisher_one_sided(freq_by_cancer, metric, "greater")) %>%
  ungroup() %>% mutate(p_adj = p.adjust(p, method="BH"))
readr::write_tsv(res_mut, "enrichment_results_mut.tsv")

res_cnv <- tibble(metric = targets_cnv) %>%
  rowwise() %>% mutate(p = fisher_one_sided(freq_by_cancer, metric, "greater")) %>%
  ungroup() %>% mutate(p_adj = p.adjust(p, method="BH"))
readr::write_tsv(res_cnv, "enrichment_results_cnv.tsv")

# ========= F. 其它通路富集（以癌种频率为单位）=========
library(msigdbr)
msig <- msigdbr(species = "Homo sapiens", category = "H") %>%
  transmute(gs_name, gene_symbol = toupper(gene_symbol)) %>%
  group_by(gs_name) %>% summarise(genes = list(unique(gene_symbol)), .groups="drop")

# 用 MAF 做“任一通路基因非同义突变”频率
#（为速度起见，这里复用 mut_freq 的计算方式：建议把 MAF 读一次缓存后循环）
# 简化起见：这里演示框架，真跑建议把每癌种 MAF 存盘后复用。
calc_path_mut_freq <- function(project, genes){
  tumor <- sub("^TCGA-", "", project)
  maf_df <- tryCatch(GDCquery_Maf(tumor = tumor, pipelines = "mutect2"),
                     error=function(e) NULL)
  if (is.null(maf_df) || nrow(maf_df)==0) return(NULL)
  maf_ns <- maf_df %>% filter(Variant_Classification %in% nonsyn) %>%
    mutate(Hugo_Symbol = toupper(Hugo_Symbol),
           barcode = norm_barcode(Tumor_Sample_Barcode))
  all_bar <- maf_ns %>% distinct(barcode) %>% pull()
  bar_hit <- maf_ns %>% filter(Hugo_Symbol %in% genes) %>% distinct(barcode) %>% pull()
  tibble(cohort=project, n=length(all_bar), hit_n=length(intersect(all_bar, bar_hit)),
         freq = ifelse(n>0, hit_n/n, NA_real_))
}

test_pathways_mut <- function(){
  out <- lapply(unique(group_df$cohort), function(prj){
    lapply(seq_len(nrow(msig)), function(i){
      gs <- msig$genes[[i]]; name <- msig$gs_name[i]
      calc_path_mut_freq(prj, gs) %>% mutate(pathway=name)
    }) %>% bind_rows()
  }) %>% bind_rows()
  out2 <- out %>% left_join(group_df, by="cohort")
  final <- out2 %>%
    group_by(pathway) %>%
    summarise(
      p = {
        x <- cur_data_all()
        # 两组 Fisher：把每癌种 freq*n 四舍五入为计数
        tbl <- x %>% group_by(group) %>%
          summarise(pos = sum(round(freq * n), na.rm=TRUE),
                    neg = sum(n - round(freq * n), na.rm=TRUE), .groups="drop")
        if (nrow(tbl)<2) NA_real_ else {
          m <- as.matrix(tbl[,c("pos","neg")]); rownames(m) <- tbl$group
          if (!all(c(low_name, high_name) %in% rownames(m))) NA_real_
          else fisher.test(m[c(low_name, high_name],), alternative="greater")$p.value
        }
      },
      .groups="drop"
    ) %>% mutate(p_adj = p.adjust(p, method="BH"))
  final
}
path_mut <- test_pathways_mut()
readr::write_tsv(path_mut, "pathway_enrichment_mut.tsv")
