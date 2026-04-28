# Manually fix some of the weirdness with the fasta (sample headers, formatting, etc)
import re
import sys


def process_fasta(input_fasta, output_fasta):
    records = []
    header = None
    seq = []

    # Parse the FASTA file
    with open(input_fasta) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append({"header": header, "seq": "".join(seq)})
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append({"header": header, "seq": "".join(seq)})

    # Standardize CF sample headers and sort by sample number
    pattern = re.compile(r'^(CF[\-_]*)(\d+)(.*)', re.IGNORECASE)  # Extract CF prefix, number, and remaining header text
    processed = []
    for rec in records:
        match = pattern.search(rec["header"])
        if match:
            num = int(match.group(2))
            rest = match.group(3)
            processed.append({"new_header": f"CF_{num:03d}{rest}", "sample_num": num, "seq": rec["seq"]})  # Reformat as CF_NNN with zero-padding to 3 digits
        else:
            processed.append({"new_header": rec["header"], "sample_num": float('inf'), "seq": rec["seq"]})

    processed.sort(key=lambda x: x["sample_num"])  # Sort by numeric sample ID, non-matching headers go last (inf)

    # Write sorted, standardized FASTA
    with open(output_fasta, 'w') as f:
        for rec in processed:
            f.write(f">{rec['new_header']}\n")
            seq_str = rec["seq"]
            for i in range(0, len(seq_str), 80):  # Wrap sequence to 80 characters per line
                f.write(f"{seq_str[i:i+80]}\n")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: fasta_sortandclean.py <input.fasta> <output.fasta>")
        sys.exit(1)
    process_fasta(sys.argv[1], sys.argv[2])
