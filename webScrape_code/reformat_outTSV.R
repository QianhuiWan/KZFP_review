

library(tidyverse)
library(readr)


kzfp_papers <- read_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/githubDirs/coh_bioLLM/tests_QW_R2/KZFP_func/outputs/kzfp_all_papers.tsv")

kzfp_otherFuntion_papers <- kzfp_papers %>% 
  # dplyr::filter(`Predicted Focus` == "Other function") %>% 
  dplyr::filter(`Other Score` >5 ) %>% 
  write_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/githubDirs/coh_bioLLM/tests_QW_R2/KZFP_func/outputs/kzfp_otherFuntion_papers.tsv")



