
# steps for processing RNA-seq counts data from TCGA database:

- step 1: get gene count matrix for all cancer types
    - ~/githubRepo/KZFP_review/tcga_KZFP/TCGACountMat_ana_S1_R2_addFilter.R

- step 2: perform DE analysis, cancer vs normal for each cancer type
    - ~/githubRepo/KZFP_review/tcga_KZFP/TCGACountMat_ana_S2_R1.R
    
- step 3: visulization
    - ~/githubRepo/KZFP_review/tcga_KZFP/TCGACountMat_ana_S3_R2_noNS.R


# steps for muation enrichment analysis

- step 1: get mutation data and CNV data
~/githubRepo/KZFP_review/tcga_KZFP/TCGA_mutationCNA_enrichemnt_R3.R
