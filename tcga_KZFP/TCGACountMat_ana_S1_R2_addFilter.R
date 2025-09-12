
# step one (S1) get raw counts and dgelist for all datasets 

.libPaths()
options(scipen=999)

library(tidyverse)
library(TCGAbiolinks)
library(edgeR)
out_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/rawData/TCGA_rnaseq"
outRes_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_tcga"

# Query for all TCGA projects with RNA-Seq data
projects <- TCGAbiolinks::getGDCprojects()
tcga_projects <- projects[grepl("TCGA", projects$project_id), "project_id"] #33

# Initialize an empty list to store projects with tumor and normal samples
# sink("TCGA_countMat_download.log")

tumor_normal_projects <- list()
tumor_normal_project_expQuery <- list()
# Loop through each project to check for both tumor and normal samples
for (project in tcga_projects) {
  query <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    # data.type = c("Aligned Reads", "Unaligned Reads"),
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = c("Primary Tumor", "Solid Tissue Normal")
    # sample.type = c("Solid Tissue Normal")
    )
  
  # If both sample types are available, add the project to the list
  # Check for available results before downloading
  results <- getResults(query)
  if (!is.null(results) && nrow(results) > 0) {
    # Download the data to the specific project directory
    # GDCdownload(query = query, directory = out_dir, files.per.chunk=4)
    
    # Add project to the list of completed downloads
    tumor_normal_projects <- c(tumor_normal_projects, project)
    tumor_normal_project_expQuery[[project]] <- query
  } else {
    message("No results found for project: ", project)
  }
}


# Display the projects with both tumor and normal samples
tumor_normal_projects #32
tumor_normal_project_expQuery #32
# sink()  # Turn off logging

saveRDS(tumor_normal_project_expQuery, 
        file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_tcga/tumor_normal_project_expQueryList_n32.rds")

tumor_normal_project_expQuery <- readRDS(paste0(outRes_dir, "/tumor_normal_project_expQueryList_n32.rds"))



# Prepare the counts files for analysis
sampInfo_dsl <- list()
countMat_dsl <- list()
raw_countMat_dsl <- list()
dgeList_dsl <- list()
for (i in 1:length(tumor_normal_project_expQuery)) {
  ## project name, i.e. name for dataset
  pro_i <- names(tumor_normal_project_expQuery)[i]
  query_i <- tumor_normal_project_expQuery[[pro_i]]
  ## sample info.
  sampInfo_res <- getResults(query_i)
  
  ## raw counts and filter
  raw_counts <- GDCprepare(query = query_i, 
                           directory = out_dir,
                           summarizedExperiment = FALSE)
  raw_counts <- raw_counts[1:60660,] # there are 60660 genes in total
  
  ## Extract log2 transformed FPKM values
  # cols <- grep("fpkm_unstranded_", colnames(raw_counts))
  # dataGBMcomplete <- log2(as.data.frame(raw_counts)[,cols]+1)
  ## Extract raw counts, not FPKM
  cols <- grep("^unstranded_TCGA", colnames(raw_counts))
  rawCounts_df <- as.data.frame(raw_counts)[,cols]
  ## Add rownames and clean the colnames
  row.names(rawCounts_df) <- paste(raw_counts$gene_name, raw_counts$gene_id, sep="_")
  colnames(rawCounts_df) <- gsub("unstranded_", "", colnames(rawCounts_df))
  
  ## create DGElist
  ### get group info.
  if (unique(sampInfo_res$cases == colnames(rawCounts_df))) {
    groups <- sampInfo_res$sample_type
  }
  dge_rawCounts <- DGEList(counts=rawCounts_df, 
                           genes=raw_counts[,1:3], #remove.zeros = TRUE,
                           group = groups)
  
  ## normalize counts with TMM method
  dge_rawCounts <- calcNormFactors(dge_rawCounts, method = "TMM")
  normalized_lcpm <- cpm(dge_rawCounts, normalized.lib.sizes = TRUE, log = TRUE)
  ############## ADD filter HERE!!!!!!!!!########################
  
  ## write sampInfo and normalized_lcpm to lists
  if (unique(sampInfo_res$cases == colnames(normalized_lcpm))) {
    sampInfo_dsl[[pro_i]] <- sampInfo_res
    countMat_dsl[[pro_i]] <- normalized_lcpm
    raw_countMat_dsl[[pro_i]] <- rawCounts_df
    dgeList_dsl[[pro_i]] <- dge_rawCounts
  }
}

# save RDSs
saveRDS(sampInfo_dsl, file = paste0(outRes_dir, "/sampInfo_dsl.rds"))
saveRDS(countMat_dsl, file = paste0(outRes_dir, "/countMat_lcpm_dsl.rds"))
saveRDS(raw_countMat_dsl, file = paste0(outRes_dir, "/raw_countMat_dsl.rds"))
saveRDS(dgeList_dsl, file = paste0(outRes_dir, "/dgeList_dsl.rds"))

# save df list to tsvs
for (name in names(sampInfo_dsl)) {
  ## save sample sheets
  dir.create(path = paste0(outRes_dir, "/sampleSheets"))
  file_path <- file.path(outRes_dir, paste0("/sampleSheets/", name, "_ss.tsv"))
  write_tsv(sampInfo_dsl[[name]], file = file_path)
  ## save count mats
  dir.create(path = paste0(outRes_dir, "/counts_lcpm"))
  file_path <- file.path(outRes_dir, paste0("/counts_lcpm/", name, "_lcpm.tsv"))
  countMat_ds_df <- as.data.frame(countMat_dsl[[name]])
  countMat_ds_df <- cbind(gene = rownames(countMat_ds_df), countMat_ds_df)
  write_tsv(countMat_ds_df, file = file_path)
  ## save raw counts
  dir.create(path = paste0(outRes_dir, "/counts_raw"))
  file_path <- file.path(outRes_dir, 
                         paste0("/counts_raw/", name, "_rawCounts.tsv"))
  raw_countMat_ds_df <- cbind(gene = rownames(raw_countMat_dsl[[name]]), 
                          raw_countMat_dsl[[name]])
  write_tsv(raw_countMat_ds_df, file = file_path)
}

