#!/usr/bin/env python3
import sys
from Bio import SeqIO

# Generated Laxtex table for proportion of each chromosome to the total genomome legth and then all other unplaced scaffolds. 
# Require fasta file to be sorted
# Accepts fasta file

def format_number(n):
    return f"{n:,}"

def main(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total_len = sum(len(rec.seq) for rec in records)

# Specific the number of chromosomes (Note: First n chromosomes must be the largest)
    chrom_records = records[:15]
    unplaced_records = records[15:]

    rows = []

    # First 15 chromosomes
    for i, rec in enumerate(chrom_records, start=1):
        seq_len = len(rec.seq)
        proportion = (seq_len / total_len) * 100
        rows.append(
            f"Chromosome {i} & {format_number(seq_len)} & {proportion:.2f}\\% \\\\"
        )

    # Combine all remaining as "Unplaced Scaffolds"
    if unplaced_records:
        unplaced_len = sum(len(r.seq) for r in unplaced_records)
        unplaced_prop = (unplaced_len / total_len) * 100
        rows.append(
            f"Unplaced Scaffolds & {format_number(unplaced_len)} & {unplaced_prop:.2f}\\% \\\\"
        )

    print("\n".join(rows))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python format_fasta_table.py <input.fasta>")
        sys.exit(1)
    main(sys.argv[1])

