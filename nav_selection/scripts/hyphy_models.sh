#!/bin/bash

set -euo pipefail

source "${CADDISRNA_ROOT:-$(pwd)}/nav_selection/config/model_config.sh"

threads="${NSLOTS:-$THREADS}"
tree_file="${IQTREE_PREFIX}.treefile"

# Replace internal stop codons with gaps so HyPhy will accept the alignment
python3 "$SCRIPT_DIR/clean_hyphy_alignment.py" "$ALIGNMENT_FILTERED" "$HYPHY_ALIGNMENT"

env CPU="$threads" hyphy busted --alignment "$HYPHY_ALIGNMENT" --tree "$tree_file" --output "$BUSTED_JSON" > "$LOG_DIR/busted.log" 2>&1
env CPU="$threads" hyphy fel    --alignment "$HYPHY_ALIGNMENT" --tree "$tree_file" --output "$FEL_JSON"    > "$LOG_DIR/fel.log"    2>&1
env CPU="$threads" hyphy meme   --alignment "$HYPHY_ALIGNMENT" --tree "$tree_file" --output "$MEME_JSON"   > "$LOG_DIR/meme.log"   2>&1
