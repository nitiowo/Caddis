#!/bin/bash --login
# usage: bash run_tblastn.sh <reference_blastdb>

module load bio/0724

query_file=$1 # query sequence

blastdb=$2 # What you are blasting against
base=$(basename $blastdb transcriptome_db)

tblastn -query $query_file \
-db $blastdb \
-out ${base}_tblastn_results.txt \
-evalue 1e-35 \
-outfmt 6