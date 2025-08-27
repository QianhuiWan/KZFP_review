
# packages needed
library(tidyverse)
library(plyranges)
library(GenomeInfoDb)

# step 1: read in the KZFP coordinate tsv file
kzfp_loc <- read_tsv("~/githubRepo/KZFP_review/input_data/kzfps_hg38_all_geneCoordinates.tsv")


# step 2: df to GR object and make sure seqnames are in UCSC format (e.g. "chr19")
kzfp_loc_GR <- kzfp_loc %>% as_granges()
seqlevelsStyle(kzfp_loc_GR) <- "UCSC"

# step 3: filter for only chr1:22 and X and Y
selected_chr <- paste0("chr", c(1:22, "X", "Y"), sep="")
kzfp_loc_GR_fil <- kzfp_loc_GR %>% 
  plyranges::filter(seqnames %in% selected_chr) %>% 
  plyranges::mutate(name = gene_name) %>% 
  plyranges::select(name)

# step 4: write into a bed file
write_bed(kzfp_loc_GR_fil, file = "~/githubRepo/KZFP_review/input_data/kzfps_hg38_all_geneCoordinates.bed")

# step 5: filter for only cancer related KZFPs and write into a bed file
# tumor suppressive n=11
tumor_suppressor <- c("ZNF23", "ZNF382", "ZNF545", "ZFP57", "ZNF471", "ZNF671", "ZBRK1",
                      "ZNF569", "ZNF496", "ZNF132", "ZNF331")
# oncogene n=14
oncogene <- c("ZNF498", "ZNF268", "ZNF224", "ZNF475", "ZNF568", "ZKSCAN3", "ZNG304",
              "ZNF10", "ZNF267", "ZNF139", "ZNF436", "ZNF689", "ZNF300", "ZNF165")

kzfp_loc_GR_fil_cancer <- kzfp_loc_GR_fil %>% 
  plyranges::filter(name %in% c(tumor_suppressor, oncogene)) %>% 
  plyranges::mutate(group = ifelse(name %in% tumor_suppressor, "repressor", "oncogene")) 

# write_bed(kzfp_loc_GR_fil_cancer, file = "~/githubRepo/KZFP_review/input_data/kzfps_hg38_cancer_geneCoordinates.bed")

# Write BED with gene names in the 4th column, and add group info in the 7th column

kzfp_loc_GR_fil_cancer_df <- kzfp_loc_GR_fil_cancer %>% 
  as.data.frame() %>% 
  dplyr::mutate(seqnames=as.character(seqnames), score = rep(0, nrow(.))) %>% 
  dplyr::select(seqnames, start, end, name, score, strand, group)

write.table(
  kzfp_loc_GR_fil_cancer_df,
  file = "~/githubRepo/KZFP_review/input_data/kzfps_hg38_cancer_geneCoordinates.bed",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)




