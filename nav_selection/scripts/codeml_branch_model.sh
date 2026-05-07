#!/bin/bash

set -euo pipefail

source "${CADDISRNA_ROOT:-$(pwd)}/nav_selection/config/model_config.sh"

# Convert FASTA to PHYLIP, label the L. flavastellus clade with #1, and write codeml tree
python3 "$SCRIPT_DIR/prepare_codeml_inputs.py" \
  "$ALIGNMENT_FILTERED" \
  "${IQTREE_PREFIX}.treefile" \
  "$FOREGROUND_TIPS_FILE" \
  "$CODEML_DIR/nav_codon.phy" \
  "$CODEML_DIR/codeml_tree.treefile"

write_ctl() {
  cat > "$CODEML_DIR/$1" <<CTL
      seqfile = nav_codon.phy
     treefile = codeml_tree.treefile
      outfile = $2
        noisy = 3
      verbose = 1
      runmode = 0
      seqtype = 1
    CodonFreq = 2
        icode = 0
        model = 2
      NSsites = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = $3
        omega = $4
    cleandata = 1
CTL
}

write_ctl codeml_alt.ctl  codeml_alt.out  0 1.5
write_ctl codeml_null.ctl codeml_null.out 1 1

cd "$CODEML_DIR"
codeml codeml_alt.ctl  > "$LOG_DIR/codeml_alt.log"  2>&1
codeml codeml_null.ctl > "$LOG_DIR/codeml_null.log" 2>&1
