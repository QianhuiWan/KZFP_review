#!/bin/bash
#SBATCH --job-name=ko_kzfpPrimate_R2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=96G
#SBATCH --time=24:00:00 # Time limit hrs:min:sec
#SBATCH --output=ko_kzfpPrimate_R2_%j.log

module load R/RStudio_R-4.4.1

# Rscript ~/githubRepo/KZFP_review/virtual_ko/scTenifoldKnk_kzfp_ko_pipe_R2test.R
Rscript ~/githubRepo/KZFP_review/virtual_ko/scTenifoldKnk_kzfp_ko_pipe_R2test_v2.R

# Rscript ~/githubRepo/KZFP_review/virtual_ko/scTenifoldKnk_kzfpPrimate_ko_pipe_R2.R


