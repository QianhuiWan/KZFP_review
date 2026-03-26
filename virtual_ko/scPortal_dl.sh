# download processed data from scPortal: https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295#study-summary

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP1162&auth_code=h0gFyD96&directory=all&context=study"  -o cfg.txt; curl -K cfg.txt && rm cfg.txt

# the other data of the above study should be in GEO with accession number GSE178341 

