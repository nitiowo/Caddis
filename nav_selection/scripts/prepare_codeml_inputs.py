# Convert codon FASTA + IQ-TREE tree into PHYLIP and codeml-ready branch-labeled tree

from pathlib import Path
import sys
from ete3 import Tree

aln_path  = Path(sys.argv[1])
tree_path = Path(sys.argv[2])
tips_path = Path(sys.argv[3])
phy_path  = Path(sys.argv[4])
tree_out  = Path(sys.argv[5])

records = []
name = None
seq = []
for line in aln_path.read_text().splitlines():
    if line.startswith(">"):
        if name is not None:
            records.append((name, "".join(seq)))
        name = line[1:].split()[0]
        seq = []
    elif line.strip():
        seq.append(line.strip())
records.append((name, "".join(seq)))

# codeml requires <=10-char sequence names, so swap full IDs for t001, t002...
name_map = {orig: f"t{i:03d}" for i, (orig, _) in enumerate(records, start=1)}

phy_path.parent.mkdir(parents=True, exist_ok=True)
seq_len = len(records[0][1])
with phy_path.open("w") as f:
    f.write(f"{len(records)} {seq_len}\n")
    for orig, s in records:
        f.write(f"{name_map[orig]:<10}  {s}\n")

# Mark the MRCA of the foreground tips with #1 for the codeml branch model
fg_tips = [t.strip() for t in tips_path.read_text().splitlines() if t.strip() and not t.startswith("#")]
tree = Tree(tree_path.read_text(), format=0)
mrca = tree.get_common_ancestor(fg_tips)
mrca.name = (mrca.name or "") + "#1"
tree_text = tree.write(format=1)

for orig in sorted(name_map, key=len, reverse=True):
    tree_text = tree_text.replace(orig, name_map[orig])
tree_out.write_text(tree_text.strip() + "\n")
