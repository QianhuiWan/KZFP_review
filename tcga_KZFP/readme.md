
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


**Notes** 
- TCGA data location: 
  - RNAseq counts: /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/TCGA_data/TCGA_rnaseq
  - SNV maf: /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/TCGA_data/TCGA_snv
  - CNV segment: /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/TCGA_data/TCGA_cnv



- File transfer command:
```
rsync -avP /Users/qwan/GDCdata qwan@apollo-acc.coh.org:/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/pub_data/TCGA_data/

```