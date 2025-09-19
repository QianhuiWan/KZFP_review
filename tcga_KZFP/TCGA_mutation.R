
# packages

library(TCGAbiolinks)
library(dplyr); library(tidyr); library(stringr); library(purrr); library(tibble); 
library(readr)
library(SummarizedExperiment)

# -------- IFN genes --------
ifn_genes <- unique(c(
  # receptors
  "IFNAR1","IFNAR2","IFNGR1","IFNGR2","IFNLR1","IL10RB",
  # JAK-STAT
  "JAK1","JAK2","TYK2","STAT1","STAT2",
  # TFs
  "IRF1","IRF3","IRF7","IRF9",
  # Tyep I IFN genes
  paste0("IFNA", 1:14), "IFNB1"
))

# nonsynonymous mutations（MAF Variant_Classification）
nonsyn_classes <- c(
  "Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del","Frame_Shift_Ins",
  "In_Frame_Del","In_Frame_Ins","Splice_Site","Translation_Start_Site","Nonstop_Mutation"
)

# -------- function：to get MAF（mutect2）--------
fetch_maf <- function(project){
  message("Downloading MAF for: ", project)
  df <- tryCatch({
    maf <- GDCquery_Maf(project, pipelines = "mutect2")
    maf %>%
      transmute(
        project = project,
        sample_id = substr(Tumor_Sample_Barcode, 1, 16),
        Hugo_Symbol,
        Variant_Classification
      )
  }, error = function(e){
    message("  ! MAF not available for ", project, " (", e$message, ")")
    tibble(project = character(), sample_id = character(), Hugo_Symbol = character(), Variant_Classification = character())
  })
  df
}

# -------- function：to get gene expression CNA（GISTIC2 cutoff）--------
fetch_cna_gene <- function(project){
  message("Downloading gene-level CNA for: ", project)
  # normal used data.type names
  data_type_candidates <- c("Gene Level Copy Number Scores", "Gene Level Copy Number")
  for (dt in data_type_candidates){
    q <- tryCatch({
      GDCquery(project = project,
               data.category = "Copy Number Variation",
               data.type     = dt)
    }, error = function(e) NULL)
    if (!is.null(q)){
      ok <- tryCatch({ GDCdownload(q); TRUE }, error = function(e) FALSE)
      if (!ok) next
      se <- tryCatch({ GDCprepare(q) }, error = function(e) NULL)
      if (is.null(se)) next

      # get marix from SummarizedExperiment 
      # 常见结构：行是基因（或含基因列名），列是样本；数值为 -2..+2
      assay_mat <- tryCatch(assay(se), error = function(e) NULL)
      if (!is.null(assay_mat)){
        cna_df <- as.data.frame(assay_mat) %>%
          tibble::rownames_to_column("Hugo_Symbol") %>%
          pivot_longer(-Hugo_Symbol, names_to = "sample_id", values_to = "CNA") %>%
          mutate(project = project)
        # get sample_id for 16 letters（TCGA-..-..-..）
        cna_df$sample_id <- substr(cna_df$sample_id, 1, 16)
        return(cna_df)
      }

      # 如果 assay 为空，试试 se@metadata 或 colData/rowData 组合（较少见）
    }
  }
  message(" ! Gene-level CNA not found for ", project)
  tibble(project = character(), sample_id = character(), Hugo_Symbol = character(), CNA = numeric())
}

# -------- 拉取所有 TCGA 项目列表 --------
all_proj <- getGDCprojects() %>%
  as.data.frame() %>%
  filter(grepl("^TCGA-", project_id)) %>%
  pull(project_id) %>%
  sort()

message("Total TCGA projects: ", length(all_proj))

# 若只想测试少量项目，解除下一行注释：
# all_proj <- c("TCGA-BRCA","TCGA-LUAD","TCGA-LAML")

# -------- 下载合并 MAF 与 CNA --------
maf_all <- map_dfr(all_proj, fetch_maf)
cna_all <- map_dfr(all_proj, fetch_cna_gene)

# -------- 生成样本状态标签 --------
# 1) TP53 非同义突变
tp53_mut <- maf_all %>%
  filter(Hugo_Symbol == "TP53", Variant_Classification %in% nonsyn_classes) %>%
  distinct(project, sample_id) %>%
  mutate(TP53_mut = 1L)

# 2) IFN 通路：非同义突变或深度 CNA（±2）
ifn_mut <- maf_all %>%
  filter(Hugo_Symbol %in% ifn_genes, Variant_Classification %in% nonsyn_classes) %>%
  distinct(project, sample_id) %>%
  mutate(IFN_mut = 1L)

# 注意：有的项目的 CNA 可能不是整数（例如阈值以外得分）；这里把 ±2 定义为“深度事件”
ifn_cna <- cna_all %>%
  filter(Hugo_Symbol %in% ifn_genes) %>%
  mutate(CNAi = suppressWarnings(as.integer(round(as.numeric(CNA))))) %>%
  filter(CNAi %in% c(-2L, 2L)) %>%
  distinct(project, sample_id) %>%
  mutate(IFN_cna_deep = 1L)

# 合并到所有样本集合（取并集）
all_samples <- bind_rows(
  maf_all %>% distinct(project, sample_id),
  cna_all %>% distinct(project, sample_id)
) %>% distinct()

res <- all_samples %>%
  left_join(tp53_mut, by = c("project","sample_id")) %>%
  left_join(ifn_mut,  by = c("project","sample_id")) %>%
  left_join(ifn_cna,  by = c("project","sample_id")) %>%
  mutate(
    TP53_mut      = ifelse(is.na(TP53_mut), 0L, 1L),
    IFN_mut       = ifelse(is.na(IFN_mut), 0L, 1L),
    IFN_cna_deep  = ifelse(is.na(IFN_cna_deep), 0L, 1L),
    IFN_any_alt   = as.integer((IFN_mut + IFN_cna_deep) > 0)
  ) %>%
  arrange(project, sample_id)

# 保存结果
readr::write_tsv(res, "TCGA_TP53_IFN_status_TCGAbiolinks.tsv")
message("Saved: TCGA_TP53_IFN_status_TCGAbiolinks.tsv  (rows = ", nrow(res), ")")








