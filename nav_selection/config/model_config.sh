#!/bin/bash

# Runtime
THREADS="${THREADS:-8}"
IQTREE_BOOTSTRAP=1000

# Paths
CADDISRNA_ROOT="${CADDISRNA_ROOT:-$(pwd)}"
NAV_DIR="$CADDISRNA_ROOT/nav_selection"
SCRIPT_DIR="$NAV_DIR/scripts"
MODEL_RAW_DIR="$NAV_DIR/outputs/model_raw"
LOG_DIR="$MODEL_RAW_DIR/logs"

ALIGNMENT_FILTERED="$NAV_DIR/outputs/alignments/nav_codon_alignment_filtered.fasta"

IQTREE_PREFIX="$MODEL_RAW_DIR/iqtree/nav_codon"
HYPHY_DIR="$MODEL_RAW_DIR/hyphy"
HYPHY_ALIGNMENT="$HYPHY_DIR/nav_codon_alignment_hyphy_ready.fasta"
BUSTED_JSON="$HYPHY_DIR/busted.json"
FEL_JSON="$HYPHY_DIR/fel.json"
MEME_JSON="$HYPHY_DIR/meme.json"

CODEML_DIR="$MODEL_RAW_DIR/codeml_branch"
FOREGROUND_TIPS_FILE="$NAV_DIR/config/foreground_tips.txt"
