# ---- Install dependencies (run once if needed) ----
# install.packages(c("tidyverse","data.table"))
# BiocManager::install(c("TCGAbiolinks","maftools"))
# install.packages("msigdbr")

library(tidyverse)
library(data.table)
library(TCGAbiolinks)
library(msigdbr)

# ========= A. Read your cancer-type grouping =========
# Expect a 2-column TSV: cohort, group
# Example rows: "TCGA-BRCA  KZFP_low", "TCGA-LUAD  KZFP_high"
group_df <- tibble(
  cancerType = c(
  "TCGA-CESC","TCGA-UCEC","TCGA-READ","TCGA-COAD","TCGA-ESCA",
  "TCGA-HNSC","TCGA-LUSC","TCGA-KIRP","TCGA-KIRC","TCGA-BLCA",
  "TCGA-STAD","TCGA-BRCA","TCGA-GBM","TCGA-KICH","TCGA-PRAD",
  "TCGA-LUAD","TCGA-PAAD","TCGA-CHOL","TCGA-LIHC","TCGA-SKCM",
  "TCGA-SARC","TCGA-THCA","TCGA-PCPG"),
  group = c(rep("KZFP_low", 15), "KZFP_high", "KZFP_noChange",
            rep("KZFP_high", 2), rep("KZFP_noChange", 4))
)

write_tsv(group_df, "~/githubRepo/KZFP_review/outputs/cancerType_groups.tsv")
  
# group_df <- readr::read_tsv("cancerType_groups.tsv") %>%
#   mutate(group = as.character(group))

# Optional: explicitly set the two group names (or keep as-is)
low_name  <- unique(group_df$group)[1]
high_name <- unique(group_df$group)[2]


# ========= B. Gene sets =========
# IFN pathway genes can be provided via 'ifn_genes.txt' (one gene per line).
# If that file doesn't exist, fall back to a default list below.
ifn_genes_file <- "ifn_genes.txt"
if (file.exists(ifn_genes_file)) {
  ifn_genes <- toupper(readr::read_lines(ifn_genes_file))
} else {
  ifn_genes <- toupper(c(
    # receptors
    "IFNAR1","IFNAR2",            # type I
    "IFNGR1","IFNGR2",            # type II
    "IFNLR1","IL10RB",            # type III
    # kinases / TFs
    "JAK1","JAK2","TYK2",
    "STAT1","STAT2","STAT3",
    "IRF9","IRF1","IRF7","IRF8",
    # negative regulators (optional, useful context)
    "SOCS1","SOCS3"
  ))
}

ifn_isg <- c(
  "ISG15","MX1","MX2","OAS1","OAS2","OAS3","OASL",
  "IFI6","IFI27","IFI44","IFI44L","IFIT1","IFIT2","IFIT3",
  "RSAD2","CMPK2","XAF1",
  "DDX58","IFIH1","IRF7","HERC5","BST2","TRIM22"
)


tp53_genes_file <- "tp53_genes.txt"
if (file.exists(tp53_genes_file)) {
  tp53_genes <- toupper(readr::read_lines(tp53_genes_file))
} else {
  tp53_genes <- toupper(c(
    "TP53","MDM2","MDM4",
    "CDKN1A","GADD45A","RRM2B","ZMAT3","DDB2","BTG2","CCNG1","PPM1D",
    "BAX","BBC3","PMAIP1","FAS","PERP","SFN","SESN1","SESN2","TP53I3","TIGAR"))
}


# Non-synonymous mutation classes (to filter MAF)
nonsyn <- c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins",
            "Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
            "Splice_Site","Translation_Start_Site")

# Normalize TCGA barcodes to first 15 chars (sample-level)
norm_barcode <- function(x) substr(x, 1, 15)



# ========= C. Mutation frequencies (from TCGAbiolinks MAF) =========
# For each TCGA project:
#   - Download Mutect2 MAF
#   - Keep non-synonymous calls
#   - Compute: fraction of samples with TP53 mutation, and with any IFN-gene mutation
calc_mut_freq_one <- function(project){
  message("Processing MAF for: ", project)
  tumor <- sub("^TCGA-", "", project)

  q <- tryCatch(GDCquery(
    project        = project,
    data.category  = "Simple Nucleotide Variation",
    data.type      = "Masked Somatic Mutation",
    workflow.type  = "Aliquot Ensemble Somatic Variant Merging and Masking"
  ), error=function(e) NULL)
  if (is.null(maf_df) || nrow(maf_df)==0) return(NULL)
  # download maf data and load in here:
  GDCdownload(q, method = "api", files.per.chunk = 20)
  maf_df <- GDCprepare(q)

  maf_ns <- maf_df %>%
    filter(Variant_Classification %in% nonsyn) %>%
    mutate(Hugo_Symbol = toupper(Hugo_Symbol),
           barcode = norm_barcode(Tumor_Sample_Barcode))
  if (nrow(maf_ns)==0) return(NULL)

  # TP53 gene carriers
  tp53_bar <- maf_ns %>% filter(Hugo_Symbol %in% tp53_genes) %>% distinct(barcode)
  # Any IFN-gene carriers
  ifn_bar  <- maf_ns %>% filter(Hugo_Symbol %in% ifn_genes) %>% distinct(barcode)

  # Denominator = all MAF samples in the project (you could switch to clinical if desired)
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

# Run for all cohorts listed in your group file
mut_list <- lapply(unique(group_df$cancerType), calc_mut_freq_one)
mut_freq <- bind_rows(mut_list)

mut_freq <- mut_freq %>% left_join(group_df, by = c("cohort"="cancerType"))

write_tsv(mut_freq, file = "~/githubRepo/KZFP_review/outputs/mut_freq_Num_sampleLevel.tsv")


# ========= D. CNV frequencies =========
# We provide a precise, gene-set–aware CNV pipeline (segment → gene overlap),
# and also keep a broad “ANY-segment” helper for sensitivity analysis.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(biomaRt)
})

# --- D0) Build GRCh38 gene coordinates via biomaRt ---
# Returns a GRanges with mcols$symbol
build_gene_coords <- function(genes, genome_version = "GRCh38"){
  genes <- unique(toupper(genes))
  # Use current Ensembl (GRCh38); pin a specific version if you need strict reproducibility
  mart <- biomaRt::useEnsembl(biomart = "genes",
                              dataset = "hsapiens_gene_ensembl",
                              version = NULL)
  tbl <- biomaRt::getBM(attributes = c("hgnc_symbol","chromosome_name",
                                       "start_position","end_position"),
                        filters = "hgnc_symbol",
                        values  = genes,
                        mart    = mart)
  # Keep only autosomes + sex chromosomes + MT
  tbl <- tbl %>% dplyr::filter(chromosome_name %in% c(as.character(1:22),"X","Y","MT"))
  GenomicRanges::GRanges(
    seqnames = tbl$chromosome_name,
    ranges   = IRanges::IRanges(start = tbl$start_position, end = tbl$end_position),
    symbol   = toupper(tbl$hgnc_symbol)
  )
}


# --- D1) Segment→gene overlap + frequency computation (core) ---
# Inputs:
#   - project: TCGA project ID (e.g., "TCGA-BRCA")
#   - cnv_seg: data.frame from GDCprepare(Copy Number Segment) with Chromosome/Start/End/Segment_Mean/Sample
#   - gene_coords: GRanges with mcols$symbol
#   - gene_set: character vector of symbols to test (e.g., "TP53" or IFN set)
# Thresholds:
#   - loss_thr (default -0.3), amp_thr (default +0.3); tighten to -0.9/+0.9 for deep events
# Output: cohort, n, loss_n, amp_n, loss_freq, amp_freq
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

  # Collapse to sample-level “any gene in the set has loss/amp”
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

# --- D2) Keep the broad “ANY-segment” helper for sensitivity analysis ---
# This treats a sample as loss/amp if ANY segment anywhere crosses the threshold.
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
    dplyr::mutate(any_loss_freq = ifelse(n>0, loss_n/n, NA_real_),
                  any_amp_freq  = ifelse(n>0, amp_n/n,  NA_real_))
}

# --- D3) Main loop: fetch CNV segments per project and compute TP53/IFN CNV frequencies ---
# 1) Build coordinates once
tp53_coords <- build_gene_coords(tp53_genes)
ifn_coords  <- build_gene_coords(ifn_genes)

# 2) Helper to fetch a project's CNV segments
get_project_cnv_segments <- function(project){
  q.cnv <- TCGAbiolinks::GDCquery(project = project,
                                  data.category = "Copy Number Variation",
                                  data.type = "Copy Number Segment")
  tryCatch({
    TCGAbiolinks::GDCdownload(q.cnv, method="api", files.per.chunk=40)
    TCGAbiolinks::GDCprepare(q.cnv)
  }, error=function(e) NULL)
}

# Compute precise (gene-set) CNV frequencies
cnv_precise_list <- lapply(unique(group_df$cancerType), function(prj){
  cnv_seg <- get_project_cnv_segments(prj)
  if (is.null(cnv_seg) || nrow(cnv_seg)==0) {
    # Return empty rows with proper column names if project has no CNV
    return(dplyr::bind_rows(
      tibble::tibble(cohort=prj, n=0, loss_n=0, amp_n=0, loss_freq=NA_real_, amp_freq=NA_real_) %>%
        dplyr::rename(tp53_loss_freq=loss_freq, tp53_amp_freq=amp_freq),
      tibble::tibble(cohort=prj, n=0, loss_n=0, amp_n=0, loss_freq=NA_real_, amp_freq=NA_real_) %>%
        dplyr::rename(ifn_loss_freq=loss_freq,  ifn_amp_freq=amp_freq)
    ))
  }
  tp53 <- get_cnv_freq_for_genes(prj, cnv_seg, tp53_coords, gene_set = tp53_genes)
  ifn  <- get_cnv_freq_for_genes(prj, cnv_seg, ifn_coords,  gene_set = ifn_genes)

  tp53 <- tp53 %>% dplyr::rename(tp53_loss_freq = loss_freq, tp53_amp_freq = amp_freq) %>% dplyr::rename(n_tp53 = n)
  ifn  <- ifn  %>% dplyr::rename(ifn_loss_freq  = loss_freq, ifn_amp_freq  = amp_freq) %>% dplyr::rename(n_ifn  = n)

  dplyr::full_join(tp53, ifn, by = "cohort")
})

cnv_precise <- dplyr::bind_rows(cnv_precise_list)

cnv_precise <- cnv_precise %>% left_join(group_df, by = c("cohort"="cancerType"))

write_tsv(cnv_precise, file = "~/githubRepo/KZFP_review/outputs/cnv_freq_Num_sampleLevel.tsv")

# Optional: also compute the broad ANY-segment frequencies as a control
# cnv_any_list <- lapply(unique(group_df$cancerType), get_cnv_seg_freq_any)
# cnv_any <- dplyr::bind_rows(cnv_any_list)

# ========= E. Merge results and run enrichment tests =========
mut_freq <- read_tsv("~/githubRepo/KZFP_review/outputs/mut_freq_Num_sampleLevel.tsv")

# Combine: your grouping + mutation frequencies + precise CNV (and optional ANY CNV)
freq_by_cancer <- group_df %>%
  dplyr::mutate(cohort = cancerType) %>% 
  dplyr::select(-group) %>% 
  # mutation-level frequencies
  left_join(mut_freq,    by="cohort") %>%  
  dplyr::select(-group) %>% 
  # TP53/IFN gene-set CNV frequencies
  left_join(cnv_precise, by="cohort") %>% 
  dplyr::mutate(group = str_replace_all(group, "KZFP_noChange", "KZFP_high"))
  # optional ANY CNV (sensitivity analysis)
  # left_join(cnv_any,     by="cohort")       

readr::write_tsv(freq_by_cancer, "~/githubRepo/KZFP_review/outputs/mut_cnv_freq_by_cancerTypes.tsv")


# One-sided Fisher on cancer-type–level frequencies:
# Convert frequencies back to counts via round(freq * n) and compare group A vs B.
fisher_one_sided <- function(tbl, col_freq, alt = "greater"){
  x <- tbl %>% 
    transmute(group, pos = round(.data[[col_freq]] * n), neg = n - round(.data[[col_freq]] * n)) %>%
    group_by(group) %>% summarise(pos=sum(pos), neg=sum(neg), .groups="drop")
  if (nrow(x) < 2) return(NA_real_)
  m <- as.matrix(x[,c("pos","neg")])
  rownames(m) <- x$group
  if (!all(c(low_name, high_name) %in% rownames(m))) return(NA_real_)
  fisher.test(m[c(low_name, high_name), ], alternative = alt)$p.value
}

# Targets for testing:
targets_mut <- c("tp53_mut_freq","ifn_mut_freq")
# Use the precise (gene-set) CNV metrics here:
targets_cnv <- c("tp53_loss_freq","tp53_amp_freq","ifn_loss_freq","ifn_amp_freq")

res_mut <- tibble(metric = targets_mut) %>%
  rowwise() %>% mutate(p = fisher_one_sided(freq_by_cancer, metric, "greater")) %>%
  ungroup() %>% mutate(p_adj = p.adjust(p, method="BH"))
readr::write_tsv(res_mut, "enrichment_results_mut.tsv")

res_cnv <- tibble(metric = targets_cnv) %>%
  rowwise() %>% mutate(p = fisher_one_sided(freq_by_cancer, metric, "greater")) %>%
  ungroup() %>% mutate(p_adj = p.adjust(p, method="BH"))
readr::write_tsv(res_cnv, "enrichment_results_cnv.tsv")

# ========= F. Optional: pathway-level mutation enrichment (Hallmark as example) =========
# For each pathway, compute per-cancer-type mutation frequency (any gene in the pathway),
# then test enrichment across your two cancer-type groups.
library(msigdbr)
msig <- msigdbr(species = "Homo sapiens", category = "H") %>%
  transmute(gs_name, gene_symbol = toupper(gene_symbol)) %>%
  group_by(gs_name) %>% summarise(genes = list(unique(gene_symbol)), .groups="drop")

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
        # Back-calculate counts per group for a 2×2 Fisher test
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
