#!/bin/bash

set -euo pipefail

source "${CADDISRNA_ROOT:-$(pwd)}/nav_selection/config/model_config.sh"

mkdir -p "$MODEL_RAW_DIR/iqtree" "$HYPHY_DIR" "$CODEML_DIR" "$LOG_DIR"

bash "$SCRIPT_DIR/iqtree_model.sh"
bash "$SCRIPT_DIR/hyphy_models.sh"
bash "$SCRIPT_DIR/codeml_branch_model.sh"
