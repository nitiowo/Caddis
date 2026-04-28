#!/bin/bash
# Usage: run from the directory containing per-sample FASTA files
# Clusters each CF*.fasta at 95% identity and writes <sample>_rep.fasta

module load bio/0724

for file in CF*.fasta; do
    base=$(basename "$file" .fasta)
    cd-hit -i "$file" -o "${base}_rep.fasta" -c 0.95 -n 5
done

module unload bio
