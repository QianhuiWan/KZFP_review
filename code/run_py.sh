#!/bin/bash
#SBATCH --job-name=KZFP_paper_scraping
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=qwan@coh.org
#SBATCH -n 2
#SBATCH -N 1
#SBATCH -p all
#SBATCH --mem=96G
#SBATCH --time=48:00:00 # Time limit hrs:min:sec
#SBATCH --output=KZFP_paper_scraping_%j.log

source /home/qwan/.bashrc
conda activate /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/miniconda3/envs/llm_env
source /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/nonCondaPkg/LLMs-from-scratch/.venv/bin/activate

# run KZFP paper scrape, R1
python /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/githubDirs/coh_bioLLM/tests_QW_R2/KZFP_func/KZFP_otherFun_webScrape_R1.py --output_dir /net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/githubDirs/coh_bioLLM/tests_QW_R2/KZFP_func/outputs

deactivate
conda deactivate

