#!/bin/bash --login
# usage: run_DROSblastn <query_file> <reference_blastdb>
#

module load bio/0724

#que_seq='/temp180/mpfrende/nvincen2/Caddis/Reference/Lim_luna/Lim-luna_candidate_sodium_channel.fasta' # query sequence

query_file=$1
blastdb=$2 # What you are blasting against
base=$(basename $blastdb _genome_db)

blastn -query $query_file \
-db $blastdb \
-out ${base}_blastn_results.txt \
-evalue 1e-35 \
-outfmt 6