# Replace internal stop codons with --- so HyPhy will accept the codon alignment

from pathlib import Path
import sys

STOPS = {"TAA", "TAG", "TGA"}

in_path = Path(sys.argv[1])
out_path = Path(sys.argv[2])

records = []
name = None
seq = []
for line in in_path.read_text().splitlines():
    if line.startswith(">"):
        if name is not None:
            records.append((name, "".join(seq)))
        name = line[1:].split()[0]
        seq = []
    elif line.strip():
        seq.append(line.strip())
records.append((name, "".join(seq)))

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as f:
    for name, s in records:
        s = s.upper()
        codons = [s[i:i + 3] for i in range(0, len(s), 3)]
        codons = ["---" if c in STOPS else c for c in codons]
        f.write(f">{name}\n{''.join(codons)}\n")
