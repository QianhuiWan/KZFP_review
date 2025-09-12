
# Step two (S2), get TopTable and filter KZFP fold change for all datasets

.libPaths()
options(scipen=999)

library(tidyverse)
library(TCGAbiolinks)
library(biomaRt)
library(edgeR)
outRes_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_tcga"

# load RDSs
sampInfo_dsl <- readRDS(paste0(outRes_dir, "/sampInfo_dsl.rds"))
countMat_dsl <- readRDS(paste0(outRes_dir, "/countMat_lcpm_dsl.rds"))
dgeList_dsl <- readRDS(paste0(outRes_dir, "/dgeList_dsl.rds"))


# get C2H2 ZFP and KZFP genes anno
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #i.e. Human genes (GRCh38.p14)
## Retrieve genes with the C2H2-type zinc finger domain (InterPro ID: IPR007087)
c2h2_zfps <- getBM(
  attributes = c("hgnc_symbol", 
                 # "ensembl_gene_id", 
                 "gene_biotype", 
                 "interpro", 
                 "description",
                 "external_gene_name"
  ),
  filters = "interpro",
  values = "IPR013087", # this is the C2H2 domain
  mart = ensembl
)

## Retrieve all genes with KRAB domain, then filter for KZFP
krab_genes <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", 
                 "interpro", "external_gene_name"),
  filters = "interpro",
  values = "IPR001909", # this is the KRAB domain
  mart = ensembl
)
## filter KZFPs from c2h2_zfps
krab_genes_fil <- krab_genes[!(krab_genes$hgnc_symbol==""), ]
kzfps <- c2h2_zfps[c2h2_zfps$hgnc_symbol%in%krab_genes_fil$hgnc_symbol,] #372


# get logFC from DE analysis
topTable_dsl <- list()
for (i in 1:length(sampInfo_dsl)) {
  ## sample info. and lcpm mat
  sampInfo_res <- sampInfo_dsl[[i]]
  normalized_lcpm <- countMat_dsl[[i]]
  dge_rawCounts <- dgeList_dsl[[i]]
  
  ## match ss with mat by sample names
  if (unique(sampInfo_res$cases == colnames(normalized_lcpm)) & 
      unique(sampInfo_res$cases == rownames(dge_rawCounts$samples))) {
    ## project name, i.e. name for dataset
    pro_i <- names(sampInfo_dsl)[i]
    sampInfo_res <- sampInfo_dsl[[pro_i]]
    normalized_lcpm <- countMat_dsl[[pro_i]]
    dge_rawCounts <- dgeList_dsl[[pro_i]]
  }else{
    print("sample sheet and data matrix not match")
  }
  
  ## DE
  if (length(unique(sampInfo_res$sample_type))==2) {
    # filter genes with lcpm>0 in at least 1 sample, lcpm<0 set to 0
    normalized_lcpm[normalized_lcpm < 0] <- 0
    normalized_lcpm_fil <- normalized_lcpm[rowSums(normalized_lcpm) > 0, ]
    
    design <- model.matrix(~0 + sample_type, data=sampInfo_res)
    colnames(design) <- 
      str_replace_all(colnames(design), 
                      c("sample_typePrimary Tumor"="sample_typeT", 
                        "sample_typeSolid Tissue Normal"="sample_typeN"))
    ### make contrast
    con <- makeContrasts(group = sample_typeT-sample_typeN, 
                         levels=design) # Tumor - Normal
    v <- voom(dge_rawCounts, design) %>%
      lmFit(design) %>%
      contrasts.fit(con) %>%
      eBayes()
    dt <- decideTests(v)
    print(summary(dt))
    ### get topTable res
    dat <- topTable(v, coef = 1, number=Inf) %>%
      # tibble::rownames_to_column(var = "gene_name") %>%
      # distinct(gene_name, .keep_all = TRUE) %>%
      dplyr::select(-t, -B, FDR = `adj.P.Val`) %>%
      mutate(Sig = ifelse(FDR < 0.05 & logFC > 2, "Up", 
                          ifelse(FDR < 0.05 & logFC < -2, "Down", "NS")))
  }else{
    print("No normal controls")
    dat <- NA
  }
  ## write dat to list for each data set
  topTable_dsl[[pro_i]] <- dat
}


# get KZFP logFC for all datasets
## remove dataset with no normal controls
topTable_dsl_fil <- topTable_dsl[!is.na(topTable_dsl)] #24
## filter KZFPs
topTable_kzfps_dsl <- list()
for (i in 1:length(topTable_dsl_fil)) {
  pro_i <- names(topTable_dsl_fil)[i]
  topTable_ds <- topTable_dsl_fil[[pro_i]]
  topTable_ds$dataset <- rep(pro_i, nrow(topTable_ds))
  topTable_ds_kzfps <- topTable_ds[topTable_ds$gene_name%in%kzfps$hgnc_symbol,]
  topTable_kzfps_dsl[[pro_i]] <- topTable_ds_kzfps # 369 
}

topTable_kzfps_ads <- topTable_kzfps_dsl %>% 
  dplyr::bind_rows(.id = "dataset") 


write_tsv(topTable_kzfps_ads, file = paste0(outRes_dir, "/topTable_kzfps_allDataSets.tsv"))

# plot KZFP logFC in heatmap for all datasets (col)

# save top table to .tsv files:
# Loop through list and write each dataframe to a .tsv file
for (name in names(topTable_dsl_fil)) {
  write.table(topTable_dsl[[name]],
              file = paste0(outRes_dir, "/topTable_files/", 
                            name, "_topTable.tsv"),
              sep = "\t",
              quote = FALSE,
              row.names = FALSE)
}



