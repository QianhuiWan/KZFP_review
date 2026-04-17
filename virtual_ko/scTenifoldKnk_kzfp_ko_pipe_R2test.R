# ============================================================
# scTenifoldKnk pipeline
# multiple genes + fixed repeat design + stable hit summary
# normal Epi vs tumor EpiT
# ============================================================
# All KO genes were evaluated using the same repeated (5 repeats) subsampling (n=3000 cells per group)
# design for normal and tumor epithelial cells, ensuring that 
# between-gene differences were not driven by inconsistent cell selection across runs.

.libPaths("/home/qwan/R/Rstudio-WithModules-WithJupyter.4.4.1.sif/libs")
suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(scTenifoldKnk)
})

# ============================================================
# 0. user settings
# ============================================================
input_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/kzfp/colon_singleCell/SCP1162"
expr_dir <- paste0(input_dir, "/expression/5fbf5a26771a5b0db8fe7a8b")
cluster_file <- paste0(input_dir, "/cluster/crc10x_tSNE_cl_global.tsv")

gene_col_in_genes_file <- 2

## all KZFPs
human_KZFPs_age_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_age_R2.tsv")
## primate-specific KZFPs
primate_specific_KZFPs_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_primate_specific_R2.tsv")

# >>> multiple KO genes here <<<
# gKO_genes <- c("ZNF141", "ZNF382", "ZNF93")
gKO_genes <- c("ZNF93")
# gKO_genes <- primate_specific_KZFPs_R2$assigned_gene

outdir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/kzfp/ko_output/scTenifoldKnk_multiKO_fixedDesign_test"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# repeated downsampling design
# n_repeats <- 5
n_repeats <- 1
n_cells_normal <- 3000
n_cells_tumor  <- 3000
base_seed <- 12355

# gene filtering
min_cells_expr <- 10

# scTenifoldKnk parameters
qc_minLSize  <- 0
qc_minCells  <- 25
nc_lambda    <- 0.1
# nc_nNet      <- 10
nc_nNet      <- 3

# nc_nCells    <- 500
nc_nCells    <- 200

nc_nComp     <- 3
nc_scaleScores <- TRUE
nc_symmetric <- FALSE
# nc_q         <- 0.05
nc_q         <- 0.9
td_K         <- 3
td_maxIter   <- 100
td_maxError  <- 1e-5
ma_nDim         <- 2

# stable hit thresholds
stable_freq_cutoff <- 0.6
top_n_per_run      <- 300
top_n_per_run      <- 50



# ============================================================
# 1. helper functions
# ============================================================
extract_barcode <- function(x) {
  sub(".*-", "", x)
}

estimate_dense_gb <- function(n_genes, n_cells, bytes_per_element = 8) {
  n_genes * n_cells * bytes_per_element / 1024^3
}

downsample_cells <- function(cell_vec, n_target, seed = 1) {
  cell_vec <- unique(cell_vec)
  if (length(cell_vec) <= n_target) return(cell_vec)
  set.seed(seed)
  sample(cell_vec, n_target, replace = FALSE)
}

filter_genes_sparse <- function(mat, min_cells_expr = 10) {
  keep <- Matrix::rowSums(mat > 0) >= min_cells_expr
  mat[keep, , drop = FALSE]
}

extract_rank_table <- function(res_obj, top_n = 300) {
  if (is.data.frame(res_obj)) {
    df <- as_tibble(res_obj)
    
    gene_col <- names(df)[tolower(names(df)) %in% c("gene", "genes", "symbol", "feature")]
    if (length(gene_col) == 0) gene_col <- names(df)[1]
    
    score_candidates <- c("score", "scores", "weight", "weights", "dr", "stat",
                          "distance", "fc", "logfc", "abs_score", "magnitude")
    score_col <- names(df)[tolower(names(df)) %in% score_candidates]
    
    if (length(score_col) == 0) {
      numeric_cols <- names(df)[sapply(df, is.numeric)]
      if (length(numeric_cols) == 0) {
        stop("Could not identify numeric score column in scTenifoldKnk result.")
      }
      score_col <- numeric_cols[1]
    } else {
      score_col <- score_col[1]
    }
    
    out <- df %>%
      transmute(
        gene = as.character(.data[[gene_col[1]]]),
        score = as.numeric(.data[[score_col]])
      ) %>%
      filter(!is.na(gene), !is.na(score)) %>%
      mutate(abs_score = abs(score)) %>%
      arrange(desc(abs_score)) %>%
      slice_head(n = top_n)
    
    return(out)
  }
  
  if (is.list(res_obj)) {
    possible_names <- names(res_obj)
    
    candidates <- c("drGenes", "rankedGenes", "result", "results", "table", "stats")
    hit <- intersect(candidates, possible_names)
    if (length(hit) > 0) {
      return(extract_rank_table(res_obj[[hit[1]]], top_n = top_n))
    }
    
    numeric_vec_names <- possible_names[sapply(res_obj, is.numeric)]
    if (length(numeric_vec_names) > 0) {
      vec <- res_obj[[numeric_vec_names[1]]]
      if (!is.null(names(vec))) {
        out <- tibble(
          gene = names(vec),
          score = as.numeric(vec)
        ) %>%
          mutate(abs_score = abs(score)) %>%
          arrange(desc(abs_score)) %>%
          slice_head(n = top_n)
        return(out)
      }
    }
  }
  
  stop("Unable to extract ranked genes from scTenifoldKnk result. Please inspect with str(res_obj).")
}

run_one_knk <- function(mat_sparse, gKO_gene, seed, group_name, repeat_id,
                        qc_minLSize, qc_minCells, nc_lambda, nc_nNet, nc_nCells,
                        nc_nComp, nc_scaleScores, nc_symmetric, nc_q,
                        td_K, td_maxIter, td_maxError, ma_nDim, top_n_per_run, gene_outdir) {
  
  message("\n==============================")
  message("Gene: ", gKO_gene, " | Group: ", group_name, " | Repeat: ", repeat_id)
  message("Sparse dim: ", nrow(mat_sparse), " genes x ", ncol(mat_sparse), " cells")
  message("Estimated dense GB: ", round(estimate_dense_gb(nrow(mat_sparse), ncol(mat_sparse)), 2))
  
  if (!(gKO_gene %in% rownames(mat_sparse))) {
    stop("gKO gene not found in matrix: ", gKO_gene)
  }
  
  set.seed(seed)
  
  mat_dense <- as.matrix(mat_sparse)
  storage.mode(mat_dense) <- "double"
  
  res <- scTenifoldKnk(
    countMatrix = mat_dense,
    gKO = gKO_gene,
    qc_minLSize = qc_minLSize,
    qc_minCells = qc_minCells,
    nc_lambda = nc_lambda,
    nc_nNet = nc_nNet,
    nc_nCells = nc_nCells,
    nc_nComp = nc_nComp,
    nc_scaleScores = nc_scaleScores,
    nc_symmetric = nc_symmetric,
    nc_q = nc_q,
    td_K = td_K,
    td_maxIter = td_maxIter,
    td_maxError = td_maxError,
    ma_nDim = ma_nDim
  )
  
  saveRDS(
    res,
    file = file.path(gene_outdir, paste0(group_name, "_repeat", repeat_id, "_", gKO_gene, "_raw.rds"))
  )
  
  rank_df <- extract_rank_table(res, top_n = top_n_per_run) %>%
    mutate(
      KO_gene = gKO_gene,
      group = group_name,
      repeat_id = repeat_id,
      rank = row_number()
    )
  
  fwrite(
    rank_df,
    file = file.path(gene_outdir, paste0(group_name, "_repeat", repeat_id, "_", gKO_gene, "_top.tsv")),
    sep = "\t"
  )
  
  rm(mat_dense, res)
  gc()
  
  return(rank_df)
}

summarize_stable_hits <- function(rank_df, n_repeats, stable_freq_cutoff = 0.6) {
  rank_df %>%
    group_by(KO_gene, group, gene) %>%
    summarise(
      n_hits = n(),
      hit_freq = n_hits / n_repeats,
      mean_rank = mean(rank, na.rm = TRUE),
      median_rank = median(rank, na.rm = TRUE),
      mean_score = mean(score, na.rm = TRUE),
      mean_abs_score = mean(abs_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(KO_gene, group, desc(hit_freq), mean_rank, desc(mean_abs_score)) %>%
    mutate(stable = hit_freq >= stable_freq_cutoff)
}

make_comparison_table <- function(stable_summary_one_gene) {
  normal_df <- stable_summary_one_gene %>%
    filter(group == "normal_Epi") %>%
    select(
      gene,
      normal_n_hits = n_hits,
      normal_hit_freq = hit_freq,
      normal_mean_rank = mean_rank,
      normal_mean_score = mean_score,
      normal_mean_abs_score = mean_abs_score,
      normal_stable = stable
    )
  
  tumor_df <- stable_summary_one_gene %>%
    filter(group == "tumor_EpiT") %>%
    select(
      gene,
      tumor_n_hits = n_hits,
      tumor_hit_freq = hit_freq,
      tumor_mean_rank = mean_rank,
      tumor_mean_score = mean_score,
      tumor_mean_abs_score = mean_abs_score,
      tumor_stable = stable
    )
  
  full_join(normal_df, tumor_df, by = "gene") %>%
    mutate(
      normal_hit_freq = dplyr::coalesce(normal_hit_freq, 0),
      tumor_hit_freq  = dplyr::coalesce(tumor_hit_freq, 0),
      normal_mean_abs_score = dplyr::coalesce(normal_mean_abs_score, 0),
      tumor_mean_abs_score  = dplyr::coalesce(tumor_mean_abs_score, 0),
      freq_diff = tumor_hit_freq - normal_hit_freq,
      abs_score_diff = tumor_mean_abs_score - normal_mean_abs_score,
      category = case_when(
        normal_stable %in% TRUE & tumor_stable %in% TRUE ~ "stable_in_both",
        normal_stable %in% TRUE & !(tumor_stable %in% TRUE) ~ "normal_specific",
        !(normal_stable %in% TRUE) & tumor_stable %in% TRUE ~ "tumor_specific",
        TRUE ~ "unstable_or_low_support"
      )
    ) %>%
    arrange(desc(abs(freq_diff)), desc(abs(abs_score_diff)))
}

# ============================================================
# 2. read input
# ============================================================
message("Reading matrix...")
mat <- Matrix::readMM(file.path(expr_dir, "matrix.mtx.gz"))

message("Reading genes...")
genes <- data.table::fread(file.path(expr_dir, "matrix.genes.tsv"), header = FALSE)

message("Reading barcodes...")
barcodes <- data.table::fread(file.path(expr_dir, "matrix.barcodes.tsv"), header = FALSE)

rownames(mat) <- make.unique(as.character(genes[[gene_col_in_genes_file]]))
colnames(mat) <- as.character(barcodes[[1]])

message("Matrix loaded: ", nrow(mat), " genes x ", ncol(mat), " cells")

message("Reading cluster annotation...")
cl <- data.table::fread(cluster_file, header = TRUE)

if (nrow(cl) >= 1 && cl$NAME[1] == "TYPE") {
  cl <- cl[-1, ]
} else if (nrow(cl) >= 2 && cl$NAME[2] == "TYPE") {
  cl <- cl[-2, ]
}

cl$barcode <- extract_barcode(cl$NAME)

# ============================================================
# 3. define groups once
# ============================================================
normal_epi <- cl %>%
  filter(ClusterMidway == "Epi")

tumor_epit <- cl %>%
  filter(ClusterMidway == "EpiT")

normal_cells_all <- intersect(unique(normal_epi$NAME), colnames(mat))
tumor_cells_all  <- intersect(unique(tumor_epit$NAME), colnames(mat))

message("Matched cells:")
message("  normal Epi = ", length(normal_cells_all))
message("  tumor EpiT = ", length(tumor_cells_all))

if (length(normal_cells_all) == 0) stop("No matched normal Epi cells.")
if (length(tumor_cells_all) == 0) stop("No matched tumor EpiT cells.")

# ============================================================
# 4. fixed repeat design
# ============================================================
message("\nBuilding fixed repeat design ...")

repeat_designs <- vector("list", n_repeats)

for (i in seq_len(n_repeats)) {
  normal_cells_use <- downsample_cells(
    normal_cells_all,
    n_target = n_cells_normal,
    seed = base_seed + i
  )
  
  tumor_cells_use <- downsample_cells(
    tumor_cells_all,
    n_target = n_cells_tumor,
    seed = base_seed + 1000 + i
  )
  
  repeat_designs[[i]] <- list(
    repeat_id = i,
    normal_cells = normal_cells_use,
    tumor_cells  = tumor_cells_use
  )
  
  fwrite(
    data.frame(barcode = normal_cells_use),
    file = file.path(outdir, paste0("repeat", i, "_normal_Epi_cells.tsv")),
    sep = "\t"
  )
  
  fwrite(
    data.frame(barcode = tumor_cells_use),
    file = file.path(outdir, paste0("repeat", i, "_tumor_EpiT_cells.tsv")),
    sep = "\t"
  )
}

# save simple summary
repeat_design_summary <- bind_rows(lapply(repeat_designs, function(x) {
  tibble(
    repeat_id = x$repeat_id,
    n_normal_cells = length(x$normal_cells),
    n_tumor_cells  = length(x$tumor_cells)
  )
}))
fwrite(
  repeat_design_summary,
  file = file.path(outdir, "fixed_repeat_design_summary.tsv"),
  sep = "\t"
)

# ============================================================
# 5. precompute matrices for each repeat design
# ============================================================
message("\nPrecomputing matrices for fixed repeat design ...")

repeat_mats <- vector("list", n_repeats)

for (i in seq_len(n_repeats)) {
  normal_cells_use <- repeat_designs[[i]]$normal_cells
  tumor_cells_use  <- repeat_designs[[i]]$tumor_cells
  
  mat_normal <- mat[, normal_cells_use, drop = FALSE]
  mat_tumor  <- mat[, tumor_cells_use,  drop = FALSE]
  
  mat_normal <- filter_genes_sparse(mat_normal, min_cells_expr = min_cells_expr)
  mat_tumor  <- filter_genes_sparse(mat_tumor,  min_cells_expr = min_cells_expr)
  
  shared_genes <- intersect(rownames(mat_normal), rownames(mat_tumor))
  mat_normal <- mat_normal[shared_genes, , drop = FALSE]
  mat_tumor  <- mat_tumor[shared_genes,  , drop = FALSE]
  
  repeat_mats[[i]] <- list(
    repeat_id = i,
    mat_normal = mat_normal,
    mat_tumor  = mat_tumor,
    shared_genes = shared_genes
  )
  
  repeat_info <- tibble(
    repeat_id = i,
    group = c("normal_Epi", "tumor_EpiT"),
    n_cells = c(ncol(mat_normal), ncol(mat_tumor)),
    n_genes = c(nrow(mat_normal), nrow(mat_tumor)),
    est_dense_gb = c(
      estimate_dense_gb(nrow(mat_normal), ncol(mat_normal)),
      estimate_dense_gb(nrow(mat_tumor), ncol(mat_tumor))
    )
  )
  
  fwrite(
    repeat_info,
    file = file.path(outdir, paste0("repeat", i, "_matrix_design.tsv")),
    sep = "\t"
  )
}

# ============================================================
# 6. run all KO genes using fixed repeat design
# ============================================================
all_gene_summaries <- list()
all_gene_comparisons <- list()

for (gKO_gene in gKO_genes) {
  message("\n#############################################")
  message("Starting KO gene: ", gKO_gene)
  message("#############################################")
  
  if (!(gKO_gene %in% rownames(mat))) {
    warning("Skipping ", gKO_gene, ": not found in expression matrix.")
    next
  }
  
  gene_outdir <- file.path(outdir, gKO_gene)
  dir.create(gene_outdir, showWarnings = FALSE, recursive = TRUE)
  
  all_rank_tables <- list()
  
  for (i in seq_len(n_repeats)) {
    message("\nRepeat ", i, " / ", n_repeats, " for ", gKO_gene)
    
    mat_normal <- repeat_mats[[i]]$mat_normal
    mat_tumor  <- repeat_mats[[i]]$mat_tumor
    shared_genes <- repeat_mats[[i]]$shared_genes
    
    if (!(gKO_gene %in% shared_genes)) {
      warning("Skipping repeat ", i, " for ", gKO_gene,
              ": gene not retained after filtering in this repeat.")
      next
    }
    
    rank_normal <- run_one_knk(
      mat_sparse = mat_normal,
      gKO_gene = gKO_gene,
      seed = base_seed + 10 * i,
      group_name = "normal_Epi",
      repeat_id = i,
      qc_minLSize = qc_minLSize,
      qc_minCells = qc_minCells,
      nc_lambda = nc_lambda,
      nc_nNet = nc_nNet,
      nc_nCells = nc_nCells,
      nc_nComp = nc_nComp,
      nc_scaleScores = nc_scaleScores,
      nc_symmetric = nc_symmetric,
      nc_q = nc_q,
      td_K = td_K,
      td_maxIter = td_maxIter,
      td_maxError = td_maxError,
      ma_nDim = ma_nDim,
      top_n_per_run = top_n_per_run,
      gene_outdir = gene_outdir
    )
    
    rank_tumor <- run_one_knk(
      mat_sparse = mat_tumor,
      gKO_gene = gKO_gene,
      seed = base_seed + 20 * i,
      group_name = "tumor_EpiT",
      repeat_id = i,
      qc_minLSize = qc_minLSize,
      qc_minCells = qc_minCells,
      nc_lambda = nc_lambda,
      nc_nNet = nc_nNet,
      nc_nCells = nc_nCells,
      nc_nComp = nc_nComp,
      nc_scaleScores = nc_scaleScores,
      nc_symmetric = nc_symmetric,
      nc_q = nc_q,
      td_K = td_K,
      td_maxIter = td_maxIter,
      td_maxError = td_maxError,
      ma_nDim = ma_nDim,
      top_n_per_run = top_n_per_run,
      gene_outdir = gene_outdir
    )
    
    all_rank_tables[[paste0("normal_", i)]] <- rank_normal
    all_rank_tables[[paste0("tumor_", i)]]  <- rank_tumor
  }
  
  if (length(all_rank_tables) == 0) {
    warning("No successful runs for ", gKO_gene)
    next
  }
  
  rank_all <- bind_rows(all_rank_tables)
  
  fwrite(
    rank_all,
    file = file.path(gene_outdir, paste0("all_repeats_ranked_top", top_n_per_run, "_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  stable_summary <- summarize_stable_hits(
    rank_df = rank_all,
    n_repeats = n_repeats,
    stable_freq_cutoff = stable_freq_cutoff
  )
  
  fwrite(
    stable_summary,
    file = file.path(gene_outdir, paste0("stable_summary_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  comparison_table <- make_comparison_table(
    stable_summary %>% filter(KO_gene == gKO_gene)
  ) %>%
    mutate(KO_gene = gKO_gene, .before = 1)
  
  fwrite(
    comparison_table,
    file = file.path(gene_outdir, paste0("normal_vs_tumor_comparison_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  fwrite(
    comparison_table %>% filter(category == "tumor_specific"),
    file = file.path(gene_outdir, paste0("tumor_specific_stable_hits_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  fwrite(
    comparison_table %>% filter(category == "normal_specific"),
    file = file.path(gene_outdir, paste0("normal_specific_stable_hits_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  fwrite(
    comparison_table %>% filter(category == "stable_in_both"),
    file = file.path(gene_outdir, paste0("stable_in_both_hits_", gKO_gene, ".tsv")),
    sep = "\t"
  )
  
  all_gene_summaries[[gKO_gene]] <- stable_summary
  all_gene_comparisons[[gKO_gene]] <- comparison_table
}

# ============================================================
# 7. combine all genes
# ============================================================
if (length(all_gene_summaries) > 0) {
  stable_summary_all <- bind_rows(all_gene_summaries)
  fwrite(
    stable_summary_all,
    file = file.path(outdir, "ALL_KO_genes_stable_summary.tsv"),
    sep = "\t"
  )
}

if (length(all_gene_comparisons) > 0) {
  comparison_all <- bind_rows(all_gene_comparisons)
  fwrite(
    comparison_all,
    file = file.path(outdir, "ALL_KO_genes_normal_vs_tumor_comparison.tsv"),
    sep = "\t"
  )
}

writeLines(capture.output(sessionInfo()),
           con = file.path(outdir, "sessionInfo.txt"))

message("\nAll done.")
message("Results saved to: ", outdir)

