#!/bin/bash
# run_trinity.sh
# Usage: ./run_trinity.sh <left_file> <right_file> <output_dir> <CPU> <MEMORY>
# This script runs Trinity. Replace variables with desired inputs if running directly.
# Combined with submit_trinity_jobs.sh this script submits multiple trinity jobs at once (one for each sample). Input variables in job script.


if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <left_file> <right_file> <output_dir> <CPU> <MEMORY>"
    exit 1
fi

LEFT_FILE="$1"
RIGHT_FILE="$2"
OUTPUT_DIR="$3"
Num_CPU="$4"
Max_Mem="$5"
trin_path="/temp180/mpfrende/nvincen2/Caddis/trinityrnaseq.v2.15.2.simg" # Path to Trinity installation

echo "Running Trinity on:"
echo "Left file: $LEFT_FILE"
echo "Right file: $RIGHT_FILE"
echo "Output directory: $OUTPUT_DIR"

# Run Trinity (adjust parameters as needed)
singularity exec -e -B /tmp "$trin_path" Trinity --seqType fq \
        --max_memory "$Max_Mem" \
        --CPU "$Num_CPU" \
        --left "${LEFT_FILE}" \
        --right "${RIGHT_FILE}" \
        --output "$OUTPUT_DIR" \
	--full_cleanup