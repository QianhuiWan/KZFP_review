# Loading single-cell data
scRNAseq <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- read.csv(scRNAseq, row.names = 1)

# loading KZFP list
## all KZFPs

## primate-specific KZFPs


# Running scTenifoldKnk
res <- scTenifoldKnk(countMatrix = scRNAseq, gKO = "G100", qc_minLSize = 0)

plotKO(res, gKO = "G100")
