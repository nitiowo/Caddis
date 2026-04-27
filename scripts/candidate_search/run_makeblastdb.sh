#!/bin/bash --login
# usage: run_DROSmakeDB <transcriptome_file>

module load bio/0724

transcriptome=$1
base=$(basename $transcriptome _L007_trinity.Trinity.fasta)

makeblastdb -in $transcriptome \
-dbtype nucl \
-out ${base}_transcriptome_db