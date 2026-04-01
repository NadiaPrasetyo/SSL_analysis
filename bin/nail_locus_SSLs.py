# script to run nail to search through SSLs HMM profiles against the locus database (fetched from pubMLST)
# info: https://github.com/TravisWheelerLab/nail
# locus: 
# SSL3: SAUR0421
# SSL7: SAUR0426
# SSL11: SAUR0420


import subprocess

subprocess.run(["nail", "search", 
    "data/6_strains/aa/SSL_aa.fasta",
    "data/pubMLST/pubMLST_2016-2026_summary_unique_6frame_proteins.fasta",
    "--tbl-out", "data/pubMLST/nail_SSL_map_results.tbl",
])
