#!/bin/bash --login
# usage: run_tblastx.sh <query_file> <transcriptome_blastdb>

module load bio/0724

# que_seq='/temp180/mpfrende/nvincen2/Caddis/Reference/Lim_luna/Lim-luna_candidate_sodium_channel.fasta'

blastdb=$2
query_file=$1
base=$(basename $blastdb _transcriptome_db)

tblastx -query $query_file \
-db $blastdb \
-out ${base}_tblastx_results.txt \
-evalue 1e-35 \
-outfmt 6
