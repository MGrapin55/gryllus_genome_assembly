#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ======================================================
# DEFAULT FASTA PATHS
# ======================================================
DEFAULT_FASTAS = [
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/RepeatMM/NUMTs/Gpenn.chr.final.fasta.masked.numt",
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/RepeatMM/NUMTs/Gfirm.chr.final.fasta.masked.numt",
    "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/GCA_046254815.1_ASM4625481v1_genomic.fna",
    "/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/full/Gbimac2.2_Xfirst_GCA_054131365.1_ASM5413136v1_genomic.fna",
]

SPECIES_NAMES = {
    "Gpenn.chr.final.fasta.masked.numt": "G.pennsylvanicus",
    "Gfirm.chr.final.fasta.masked.numt": "G.firmus",
    "GCA_046254815.1_ASM4625481v1_genomic.fna": "G.assimilis",
    "Gbimac2.2_Xfirst_GCA_054131365.1_ASM5413136v1_genomic.fna": "G.bimaculatus",
}

COLORS = {
    "G.pennsylvanicus": "#618B4A",
    "G.firmus": "#AFBC88",
    "G.assimilis": "#7AA095",
    "G.bimaculatus": "#49306B",
}

# ======================================================
# HELPERS
# ======================================================
def format_number(n):
    return f"{n:,}"

def count_lowercase(seq):
    return sum(1 for b in seq if b.islower())

def load_softmask_stats(path):
    """
    Returns:
      total_len,
      chrom_soft (len 15),
      chrom_total (len 15),
      unplaced_soft,
      unplaced_total
    """
    records = list(SeqIO.parse(path, "fasta"))

    chrom = records[:15]
    unplaced = records[15:]

    chrom_soft = []
    chrom_total = []

    for rec in chrom:
        soft = count_lowercase(rec.seq)
        total = len(rec.seq)
        chrom_soft.append(soft)
        chrom_total.append(total)

    unplaced_soft = sum(count_lowercase(r.seq) for r in unplaced)
    unplaced_total = sum(len(r.seq) for r in unplaced)

    total_len = sum(chrom_total) + unplaced_total

    return (
        total_len,
        chrom_soft,
        chrom_total,
        unplaced_soft,
        unplaced_total
    )

# ======================================================
# PLOTTING
# ======================================================
def make_plots(df):
    # ---- SOFTMASK BP ----
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df,
        x="Region",
        y="Softmask_BP",
        hue="Species",
        dodge=True,
        linewidth=1,
        edgecolor="black",
        palette=COLORS
    )
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Chromosome")
    plt.ylabel("Softmasked base pairs")
    plt.title("Softmasked Base Pairs by Chromosome")
    plt.tight_layout()
    plt.savefig("softmask_bp.png", dpi=300)
    plt.close()

    # ---- SOFTMASK PROPORTION ----
    plt.figure(figsize=(14, 6))
    sns.barplot(
        data=df,
        x="Region",
        y="Softmask_Prop",
        hue="Species",
        dodge=True,
        linewidth=1,
        edgecolor="black",
        palette=COLORS
    )
    plt.xticks(rotation=45, ha="right")
    plt.xlabel("Chromosome")
    plt.ylabel("Softmasked proportion (%)")
    plt.title("Softmasked Proportion by Chromosome")
    plt.tight_layout()
    plt.savefig("softmask_proportion.png", dpi=300)
    plt.close()

# ======================================================
# MAIN
# ======================================================
def main(fastas):
    summaries = []
    plot_rows = []

    # ---- LOAD DATA ----
    for fa in fastas:
        (
            total_len,
            chrom_soft,
            chrom_total,
            unplaced_soft,
            unplaced_total
        ) = load_softmask_stats(fa)

        species = SPECIES_NAMES.get(os.path.basename(fa), os.path.basename(fa))

        summaries.append(
            (total_len, chrom_soft, chrom_total, unplaced_soft, unplaced_total)
        )

        # X + chr1–14
        for i in range(15):
            region = "X" if i == 0 else str(i)
            soft = chrom_soft[i]
            total = chrom_total[i]

            plot_rows.append({
                "Region": region,
                "Softmask_BP": soft,
                "Softmask_Prop": (soft / total) * 100 if total > 0 else 0,
                "Species": species
            })

        # Unplaced
        plot_rows.append({
            "Region": "Unplaced",
            "Softmask_BP": unplaced_soft,
            "Softmask_Prop": (unplaced_soft / unplaced_total) * 100 if unplaced_total > 0 else 0,
            "Species": species
        })

    # ---- LATEX TABLE ----
    for i in range(16):
        row = []

        if i == 0:
            row.append("X Chromosome")
        elif i < 15:
            row.append(f"Chromosome {i}")
        else:
            row.append("Unplaced Scaffolds")

        for (
            total_len,
            chrom_soft,
            chrom_total,
            unplaced_soft,
            unplaced_total
        ) in summaries:

            if i < 15:
                soft = chrom_soft[i]
                total = chrom_total[i]
            else:
                soft = unplaced_soft
                total = unplaced_total

            prop = (soft / total) * 100 if total > 0 else 0
            row.append(format_number(soft))
            row.append(f"{prop:.2f}\\%")

        print(" & ".join(row) + " \\\\")

    # ---- PLOTS ----
    df = pd.DataFrame(plot_rows)
    chrom_order = ["X"] + [str(i) for i in range(1, 15)] + ["Unplaced"]
    df["Region"] = pd.Categorical(df["Region"], categories=chrom_order, ordered=True)
    
    # =======================
    # WRITE TSV  ### NEW ###
    # =======================
    out_tsv = "softmasked_metrics_summary.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Wrote summary table to: {out_tsv}")

    make_plots(df)

# ======================================================
# ENTRY POINT
# ======================================================
if __name__ == "__main__":

    if len(sys.argv) > 1:
        fasta_files = sys.argv[1:]
    else:
        fasta_files = DEFAULT_FASTAS

    for f in fasta_files:
        if not os.path.exists(f):
            sys.exit(f"ERROR: FASTA not found: {f}")

    main(fasta_files)
