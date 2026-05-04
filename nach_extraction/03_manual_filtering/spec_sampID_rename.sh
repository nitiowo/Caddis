#!/bin/bash
# Usage: spec_sampID_rename.sh <species_code> <sampleID> <input_fasta>
# Prepends species code and sample ID to each FASTA header

species=$1
sampleID=$2
input=$3

awk -v sp="$species" -v sid="$sampleID" '
  /^>/ {
    print ">" sp "_" sid "_" substr($0, 2)
    next
  }
  { print }
' "$input" > "${sampleID}_hits_renamed.fasta"
