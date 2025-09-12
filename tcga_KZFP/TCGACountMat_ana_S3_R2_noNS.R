
# Step two (S3), plot KZFP fold change for all datasets

.libPaths()
options(scipen=999)

library(tidyverse)
library(TCGAbiolinks)
library(pheatmap)
library(ComplexHeatmap)
library(cluster)
library(seriation)
library(stats)
library(RColorBrewer)
library(ggplot2)
library(ggplotify)

outRes_dir <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/projects/PMM/RNA_tcga"

# load tsv
topTable_kzfps_ads <- read_tsv(paste0(outRes_dir, "/topTable_kzfps_allDataSets.tsv"))


# plot KZFP logFC in heatmap for all datasets (as columns)

topTable_kzfps_ads_mat <- 
  topTable_kzfps_ads %>% 
  # dplyr::filter(!Sig=="NS") %>%
  # dplyr::filter(abs(logFC)>=1) %>%
  dplyr::filter(abs(logFC)>=1 & FDR < 0.05) %>%
  pivot_wider(id_cols = gene_name, names_from = dataset, 
              values_from = logFC, values_fill = 0) %>% 
  tibble::column_to_rownames(var="gene_name") %>% 
  as.matrix()

# check distribution and define pan-cancer KZFPs
topTable_kzfps_ads_mat_df <- as.data.frame(topTable_kzfps_ads_mat) %>%
  mutate(
    pos_count = rowSums(across(everything(), ~ .x > 0), na.rm = TRUE),
    neg_count = rowSums(across(everything(), ~ .x < 0), na.rm = TRUE),
    all_de_count = pos_count + neg_count
  )

de_KZFP_distribution <- topTable_kzfps_ads_mat_df %>% 
  tibble::rownames_to_column(var = "KZFPs") %>% 
  dplyr::select(KZFPs, 
                up_DE = pos_count, 
                down_DE = neg_count, all_DE = all_de_count) %>% 
  pivot_longer(cols = -c(KZFPs), 
               names_to = "DE", values_to = "Num_cancerTypes") %>% 
  ggplot(aes(x = Num_cancerTypes, fill = DE))+
  geom_bar(position = position_dodge(preserve = 'single'))+
  # geom_density()+
  facet_wrap(~DE, ncol = 1, labeller = as_labeller(
    c(`all_DE` = "All DE KZFPs", `up_DE` = "Upregulated KZFPs", 
      `down_DE` = "Downregulated KZFPs")))+
  ggpubr::theme_pubr(base_size = 12)+
  labs(x = "#Cancer types", y = "#DE KZFPs")+
  theme(legend.position = "none")
  
ggsave(filename = "de_KZFP_distribution.pdf", 
       plot = de_KZFP_distribution, device = "pdf", 
       path = "/home/qwan/githubRepo/KZFP_review/outputs/", 
       width = 12, height = 10, units = "cm")


# pan-cancer KZFP mat:
keep <- (rowSums(topTable_kzfps_ads_mat > 0, na.rm = TRUE) >= 3) | (rowSums(topTable_kzfps_ads_mat < 0, na.rm = TRUE) >= 3)
panCancer_kzfps_ads_mat <- topTable_kzfps_ads_mat[keep, , drop = FALSE]


keep <- (rowSums(panCancer_kzfps_ads_mat > 0, na.rm = TRUE) >= 3)
panCancer_up_kzfps_ads_mat <- panCancer_kzfps_ads_mat[keep, , drop = FALSE]

# keep <- (rowSums(panCancer_kzfps_ads_mat < 0, na.rm = TRUE) >= 3)
keep <- (rowSums(panCancer_kzfps_ads_mat < 0, na.rm = TRUE) >= 5)
panCancer_down_kzfps_ads_mat <- panCancer_kzfps_ads_mat[keep, , drop = FALSE]

# cancer-type specific KZFP mat
keep <- (rowSums(topTable_kzfps_ads_mat > 0, na.rm = TRUE) ==1 )
typeSpecific_up_kzfps_ads_mat <- topTable_kzfps_ads_mat[keep, , drop = FALSE]

keep <- (rowSums(topTable_kzfps_ads_mat < 0, na.rm = TRUE) == 1)
typeSpecific_down_kzfps_ads_mat <- topTable_kzfps_ads_mat[keep, , drop = FALSE]


## annotations for columns and plot
### Get TCGA project metadata
projects <- TCGAbiolinks::getGDCprojects()
### Filter for TCGA projects and select relevant columns, get cancer type names
tcga_metadata <- projects[grepl("TCGA", projects$project_id), 
                          c("project_id", "name")] %>% 
  # Select relevant columns
  dplyr::select(Project = project_id, Cancer_type = name) %>%
  # Filter for projects in tcga_projects
  dplyr::filter(Project %in% colnames(topTable_kzfps_ads_mat)) %>%
  # Reorder to match tcga_projects
  arrange(match(colnames(topTable_kzfps_ads_mat), Project)) 

### add original tisue type to cancer type
tcga_metadata$TissueType <- 
  # c("Endocrine", "Digestive", "Digestive", "Soft Tissue", 
  #   "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory", 
  #   "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive", 
  #   "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain", 
  #   "Endocrine", "Skin", "Reproductive")
  # c("Endocrine", "Lymphatic", "Digestive", "Digestive", "Soft Tissue", 
  #   "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory", 
  #   "Liver", "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive", 
  #   "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain", 
  #   "Endocrine", "Skin", "Reproductive")
  c("Endocrine", "Digestive", "Digestive", "Soft Tissue",
    "Reproductive", "Digestive", "Urinary", "Head and Neck", "Respiratory",
    "Liver", "Respiratory", "Urinary", "Urinary", "Urinary", "Digestive",
    "Reproductive", "Digestive", "Reproductive", "Digestive", "Brain",
    "Endocrine", "Skin", "Reproductive")

tcga_metadata$TissueType <- factor(tcga_metadata$TissueType, 
                                   levels = unique(tcga_metadata$TissueType))



# Column annotations (e.g., Sample Type and Condition)
col_annotation <- HeatmapAnnotation(
  TissueType = tcga_metadata$TissueType,
  col = list(
    TissueType = setNames(
      colorRampPalette(brewer.pal(8, "Set3"))(11)[-2], 
      levels(tcga_metadata$TissueType))
  )
)

# Create the heatmap with column annotations
library(circlize)
color_scale <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))

## test different clustering method
# register_DendSer()
# get_seriation_method("dist", "DendSer")
methods <- list_seriation_methods()
# for (method in methods$dist[c(1,5:8, 11:length(methods$dist))]) {
for (method in methods$dist[c(5:8,11:length(methods$dist))]) {
# for (method in methods$dist[c(1)]) {
  o1 = seriate(dist(panCancer_up_kzfps_ads_mat), method = method)
  o2 = seriate(dist(t(panCancer_up_kzfps_ads_mat)), method = method)
  
  p <- Heatmap(
    panCancer_up_kzfps_ads_mat,
    top_annotation = col_annotation,
    name = "LogFC", 
    # clustering_distance_columns = "spearman",
    # clustering_method_columns = "complete",
    # clustering_distance_rows = "maximum",
    # clustering_method_rows= "ward",
    row_order = get_order(o1),
    column_order = get_order(o2),
    # cluster_rows = as.dendrogram(o1[[1]]), 
    # cluster_columns = as.dendrogram(o2[[1]]),
    column_names_centered = TRUE, 
    show_row_names = TRUE,
    col = color_scale
    # cluster_rows = diana(topTable_kzfps_ads_mat),
    # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
  )
  
  # Convert to ggplot object
  ggplot_heatmap <- as.ggplot(p)
  ggsave(paste0("/home/qwan/githubRepo/KZFP_review/tcga_KZFP/", 
                "pdf_test/kzfps_FC_TCGAds_panCancer_up_", method, ".pdf"), 
         plot = ggplot_heatmap,
         width = 20, height = 15, units = "cm", limitsize = FALSE)
}

# the heatmap for up-regulated KZFPs in pan-cancer
o1 = seriate(dist(panCancer_up_kzfps_ads_mat), method = "MDS")
o2 = seriate(dist(t(panCancer_up_kzfps_ads_mat)), method = "MDS")
# o1 = seriate(dist(topTable_kzfps_ads_mat), method = "R2E")
# o2 = seriate(dist(t(topTable_kzfps_ads_mat)), method = "R2E")
p <- Heatmap(
  panCancer_up_kzfps_ads_mat,
  # top_annotation = col_annotation,
  name = "LogFC", 
  # clustering_distance_columns = "spearman",
  # clustering_method_columns = "complete",
  # clustering_distance_rows = "maximum",
  # clustering_method_rows= "complete",
  row_order = get_order(o1),
  column_order = get_order(o2),
  # cluster_rows = as.dendrogram(o1[[1]]),
  # cluster_columns = rev(as.dendrogram(o2[[1]])),
  column_names_centered = TRUE, 
  show_row_names = TRUE,
  col = color_scale
  # cluster_rows = diana(topTable_kzfps_ads_mat),
  # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
)

# Convert to ggplot object
ggplot_heatmap <- as.ggplot(p)
ggsave(paste0("/home/qwan/githubRepo/KZFP_review/outputs/", 
              "/panCancer_up_kzfps_FCabove1_FDR05_TCGAds_MDS.pdf"), 
       plot = ggplot_heatmap,
       width = 20, height = 20, units = "cm", limitsize = FALSE)



# for down regulated pan-cancer KZFPs
## test different clustering method
for (method in methods$dist[c(1, 5:8, 11:length(methods$dist))]) {
# for (method in methods$dist[c(5:8,11:length(methods$dist))]) {
  # for (method in methods$dist[c(1)]) {
  o1 = seriate(dist(panCancer_down_kzfps_ads_mat), method = method)
  o2 = seriate(dist(t(panCancer_down_kzfps_ads_mat)), method = method)
  
  p <- Heatmap(
    panCancer_down_kzfps_ads_mat,
    top_annotation = col_annotation,
    name = "LogFC", 
    # clustering_distance_columns = "spearman",
    # clustering_method_columns = "complete",
    # clustering_distance_rows = "maximum",
    # clustering_method_rows= "ward",
    row_order = get_order(o1),
    column_order = get_order(o2),
    # cluster_rows = as.dendrogram(o1[[1]]), 
    # cluster_columns = as.dendrogram(o2[[1]]),
    column_names_centered = TRUE, 
    show_row_names = TRUE,
    col = color_scale
    # cluster_rows = diana(topTable_kzfps_ads_mat),
    # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
  )
  
  # Convert to ggplot object
  ggplot_heatmap <- as.ggplot(p)
  ggsave(paste0("/home/qwan/githubRepo/KZFP_review/tcga_KZFP/", 
                "pdf_test/kzfps_FC_TCGAds_panCancer_down_", method, "_v2.pdf"), 
         plot = ggplot_heatmap, #scale = 0.8,
         width = 20, height = 30, units = "cm", limitsize = FALSE)
}



# the heatmap for up-regulated KZFPs in pan-cancer
o1 = seriate(dist(panCancer_down_kzfps_ads_mat), method = "MDS")
o2 = seriate(dist(t(panCancer_down_kzfps_ads_mat)), method = "MDS")
# o1 = seriate(dist(topTable_kzfps_ads_mat), method = "R2E")
# o2 = seriate(dist(t(topTable_kzfps_ads_mat)), method = "R2E")
p <- Heatmap(
  panCancer_down_kzfps_ads_mat,
  # top_annotation = col_annotation,
  name = "LogFC", 
  # clustering_distance_columns = "spearman",
  # clustering_method_columns = "complete",
  # clustering_distance_rows = "maximum",
  # clustering_method_rows= "complete",
  row_order = get_order(o1),
  column_order = rev(get_order(o2)),
  # cluster_rows = as.dendrogram(o1[[1]]),
  # cluster_columns = rev(as.dendrogram(o2[[1]])),
  column_names_centered = TRUE, 
  show_row_names = TRUE,
  col = color_scale
  # cluster_rows = diana(topTable_kzfps_ads_mat),
  # cluster_columns = agnes(t(topTable_kzfps_ads_mat))
)

# Convert to ggplot object
ggplot_heatmap <- as.ggplot(p)
ggsave(paste0("/home/qwan/githubRepo/KZFP_review/outputs/", 
              "/panCancer_down_kzfps_FCabove1_FDR05_TCGAds_MDS.pdf"), 
       plot = ggplot_heatmap,
       width = 20, height = 30, units = "cm", limitsize = FALSE)


