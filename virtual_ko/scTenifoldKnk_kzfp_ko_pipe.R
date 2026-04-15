
library(tidyverse)
library(scTenifoldKnk)
library(Matrix)

# Loading single-cell data
input_dir = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/kzfp/colon_singleCell/SCP1162/expression/5fbf5a26771a5b0db8fe7a8b"
mat <- readMM(paste0(input_dir, "/matrix.mtx.gz"))
genes <- read.delim(paste0(input_dir, "/matrix.genes.tsv"), header = FALSE)
cells <- read.delim(paste0(input_dir, "/matrix.barcodes.tsv"), header = FALSE)

rownames(mat) <- genes$V2   # or V2 depending on file
colnames(mat) <- cells$V1

# load cell types
cell_types_info <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/kzfp/colon_singleCell/SCP1162/cluster/crc10x_tSNE_cl_global.tsv")
cell_types_info <- cell_types_info[-1, ]



# dgCMatrix sparse matrix to normal matrix
scRNAseq <- as.matrix(mat) # too big

# loading KZFP gene lists
## all KZFPs
human_KZFPs_age_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_age_R2.tsv")

## primate-specific KZFPs
primate_specific_KZFPs_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_primate_specific_R2.tsv")

# Running scTenifoldKnk
res <- scTenifoldKnk(countMatrix = scRNAseq, gKO = "G10", qc_minLSize = 0)

plotKO(res, gKO = "G10")



