#!/usr/bin/env python3
import sys
import os
from Bio import SeqIO

# =========================
# DEFAULT FASTA FILES HERE
# =========================
DEFAULT_FASTAS = [
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/RepeatMM/NUMTs/Gpenn.chr.final.fasta.masked.numt",
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/RepeatMM/NUMTs/Gfirm.chr.final.fasta.masked.numt",
    "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/GCA_046254815.1_ASM4625481v1_genomic.fna",
    "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/GCA_965638035.1_iqGryBima1.hap1.1_genomic.fna.NoMito",
]

def format_number(n):
    return f"{n:,}"

def summarize_fasta(fasta_path, n_chrom=15):
    """
    Returns:
      total_len: total assembly length
      chrom_lens: list of chromosome lengths (X + chr1–14)
      unplaced_len: combined length of all remaining sequences
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total_len = sum(len(r.seq) for r in records)

    chrom_records = records[:n_chrom]
    unplaced_records = records[n_chrom:]

    chrom_lens = [len(r.seq) for r in chrom_records]
    unplaced_len = sum(len(r.seq) for r in unplaced_records)

    return total_len, chrom_lens, unplaced_len


def main(fasta_files):
    summaries = [summarize_fasta(f) for f in fasta_files]

    n_chrom = len(summaries[0][1])  # number of chromosomes (X + autosomes)

    # Build rows
    for i in range(n_chrom):
        if i == 0:
            label = "X Chromosome"
        else:
            label = f"Chromosome {i}"

        row = [label]

        for total_len, chrom_lens, _ in summaries:
            seq_len = chrom_lens[i]
            prop = (seq_len / total_len) * 100
            row.append(f"{format_number(seq_len)} & {prop:.2f}\\%")

        print(" & ".join(row) + " \\\\")

    # Unplaced scaffolds row
    row = ["Unplaced Scaffolds"]
    for total_len, _, unplaced_len in summaries:
        prop = (unplaced_len / total_len) * 100
        row.append(f"{format_number(unplaced_len)} & {prop:.2f}\\%")

    print(" & ".join(row) + " \\\\")

if __name__ == "__main__":
#    if len(sys.argv) < 2:
#        print("Usage: python format_fasta_table_multi.py <assembly1.fasta> <assembly2.fasta> [...]")
#        sys.exit(1)

#    main(sys.argv[1:])
# Use command-line FASTAs if provided, otherwise defaults
    if len(sys.argv) > 1:
        fasta_files = sys.argv[1:]
    else:
        fasta_files = DEFAULT_FASTAS

    # Safety check (NOW fasta_files exists)
    for f in fasta_files:
        if not os.path.exists(f):
            sys.exit(f"ERROR: FASTA not found: {f}")

    main(fasta_files)
