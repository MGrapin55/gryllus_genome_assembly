#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ======================================================
# DEFAULT FASTA PATHS (UNCHANGED – EDIT IF NEEDED)
# ======================================================
DEFAULT_FASTAS = [
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/RepeatMM/NUMTs/Gpenn.chr.final.fasta.masked.numt",
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/RepeatMM/NUMTs/Gfirm.chr.final.fasta.masked.numt",
    "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/GCA_046254815.1_ASM4625481v1_genomic.fna",
    "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/GCA_965638035.1_iqGryBima1.hap1.1_genomic.fna.NoMito",
]
# ======================================================
# SPECIES NAME + COLOR MAP
# ======================================================
SPECIES_NAMES = {
    "Gpenn.chr.final.fasta.masked.numt": "G.pennsylvanicus",
    "Gfirm.chr.final.fasta.masked.numt": "G.firmus",
    "GCA_046254815.1_ASM4625481v1_genomic.fna": "G.assimilis",
    "GCA_965638035.1_iqGryBima1.hap1.1_genomic.fna.NoMito": "G.bimaculatus",
}

COLORS = {
    "G.pennsylvanicus": "#618B4A",
    "G.firmus": "#AFBC88",
    "G.assimilis": "#7AA095",
    "G.bimaculatus": "#49306B",
}

# ======================================================
# FASTA SUMMARY
# ======================================================
def summarize_fasta(fasta_path, n_chrom=15):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total_len = sum(len(r.seq) for r in records)

    chrom_records = records[:n_chrom]
    unplaced_records = records[n_chrom:]

    chrom_lens = [len(r.seq) for r in chrom_records]
    unplaced_len = sum(len(r.seq) for r in unplaced_records)

    return total_len, chrom_lens, unplaced_len

# ======================================================
# PLOTTING
# ======================================================
def make_plots(df):
    # --- BP plot ---
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df,
        x="Region",
        y="BP",
        hue="Species",
        dodge=True,
        linewidth=1,
        edgecolor="black",
        palette=COLORS
    )

    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Chromosome")
    plt.ylabel("Base pairs")
    plt.title("Chromosome Length by Species")
    plt.tight_layout()
    plt.savefig("chromosome_lengths_bp.png", dpi=300)
    plt.close()

    # --- Proportion plot ---
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df,
        x="Region",
        y="Proportion",
        hue="Species",
        dodge=True,
        linewidth=1,
        edgecolor="black",
        palette=COLORS
    )

    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Chromosome")
    plt.ylabel("Genome proportion (%)")
    plt.title("Chromosome Proportion by Species")
    plt.tight_layout()
    plt.savefig("chromosome_proportion_percent.png", dpi=300)
    plt.close()

# ======================================================
# MAIN
# ======================================================
def main(fasta_files):
    rows = []

    for fasta in fasta_files:
        total_len, chrom_lens, unplaced_len = summarize_fasta(fasta)
        species = SPECIES_NAMES.get(os.path.basename(fasta), os.path.basename(fasta))

        # Chromosomes (X + 1–14)
        for i, bp in enumerate(chrom_lens):
            region = "X" if i == 0 else str(i)
            rows.append({
                "Region": region,
                "BP": bp,
                "Proportion": bp / total_len * 100,
                "Species": species
            })

        # Unplaced scaffolds
        rows.append({
            "Region": "Unplaced",
            "BP": unplaced_len,
            "Proportion": unplaced_len / total_len * 100,
            "Species": species
        })

    df = pd.DataFrame(rows)

    # Preserve chromosome order
    chrom_order = ["X"] + [str(i) for i in range(1, 15)] + ["Unplaced"]
    df["Region"] = pd.Categorical(df["Region"], categories=chrom_order, ordered=True)

    make_plots(df)

# ======================================================
# ENTRY POINT
# ======================================================
if __name__ == "__main__":

    # Use CLI FASTAs if provided, otherwise defaults
    if len(sys.argv) > 1:
        fasta_files = sys.argv[1:]
    else:
        fasta_files = DEFAULT_FASTAS

    # Safety check
    for f in fasta_files:
        if not os.path.exists(f):
            sys.exit(f"ERROR: FASTA not found: {f}")

    main(fasta_files)

