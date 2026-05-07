#!/bin/bash

set -euo pipefail

source "${CADDISRNA_ROOT:-$(pwd)}/nav_selection/config/model_config.sh"

iqtree2 \
  -s "$ALIGNMENT_FILTERED" \
  -st CODON \
  -m MFP \
  -bb "$IQTREE_BOOTSTRAP" \
  -nt "${NSLOTS:-$THREADS}" \
  -pre "$IQTREE_PREFIX" \
  > "$LOG_DIR/iqtree.log" 2>&1
