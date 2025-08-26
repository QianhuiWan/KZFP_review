#!/bin/bash  
# This script sets up the environment for KZFP review analysis

## create env for kzfp review paper analysis
conda create -n kzfp_review python=3.11 -y

conda activate kzfp_review

## install py packages needed
conda install -c conda-forge numpy pandas matplotlib seaborn jupyterlab ipykernel -y
### pip is only for install python package that are not in or not 
### frequently updated in conda
pip install plotly

## setup jupyter kernel (location path should be: /Users/qwan/.local/share/jupyter/kernels/kzfp_review)
python -m ipykernel install --user --name kzfp_review --display-name "Python(kzfp_review)"

## start jupyter notebook with this command:
jupyter lab


