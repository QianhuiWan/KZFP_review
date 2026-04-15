
library(tidyverse)
library(scTenifoldKnk)

# Loading single-cell data
scRNAseq <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- read.csv(scRNAseq, row.names = 1)

# loading KZFP list
## all KZFPs
human_KZFPs_age_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_age_R2.tsv")

## primate-specific KZFPs
primate_specific_KZFPs_R2 <- read_tsv(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/KZFPs/KZFP_evolutional_primate_specific_R2.tsv")

# Running scTenifoldKnk
res <- scTenifoldKnk(countMatrix = scRNAseq, gKO = "G100", qc_minLSize = 0)

# plotKO(res, gKO = "G10")

# plotKO_fix function
plotKO_fix <- function(res, gKO) {
  
  gList <- unique(c(
    gKO,
    res$diffRegulation$gene[
      res$diffRegulation$distance > 1e-10 &
        res$diffRegulation$p.adj < 0.05
    ]
  ))
  
  if (length(gList) == 0) gList <- gKO
  
  # core fix
  sCluster <- as.matrix(res$tensorNetworks$WT[gList, gList, drop = FALSE])
  
  rownames(sCluster) <- gList
  colnames(sCluster) <- gList
  
  # only 1 gene, plot a dot
  if (length(gList) == 1) {
    plot(1, 1,
         pch = 19,
         cex = 3,
         axes = FALSE,
         xlab = "",
         ylab = "",
         main = paste("KO:", gKO))
    text(1, 1, labels = gKO, pos = 3)
    return()
  }
  
  # simple network
  library(igraph)
  
  sCluster[abs(sCluster) < quantile(abs(sCluster), 0.99)] <- 0
  diag(sCluster) <- 0
  
  net <- graph_from_adjacency_matrix(sCluster, mode = "directed", weighted = TRUE)
  
  plot(net,
       vertex.size = 10,
       vertex.label.cex = 0.7,
       edge.arrow.size = 0.3,
       main = paste("KO:", gKO))
}

plotKO_fix(res, gKO = "G100")


