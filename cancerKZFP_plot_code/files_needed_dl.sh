#!/bin/bash 
cd ~/githubRepo/KZFP_review/input_data

curl http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz -o hg38_cytoBand.txt.gz

pigz -dk hg38_cytoBand.txt.gz


