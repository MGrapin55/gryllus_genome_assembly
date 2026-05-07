#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# -------------------------
# INPUTS
# -------------------------

# Example usage:
# files = ["sp1.annotations.txt", "sp2.annotations.txt"]
# species = ["Species1", "Species2"]

files = [
    "/mnt/nrdstor/moorelab/mgrapin2/GryllusGenomeBackup/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/funannotate/funannotate/annotate_results/Gryllus_pennsylvanicus.annotations.txt",
    "/mnt/nrdstor/moorelab/mgrapin2/GryllusGenomeBackup/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/Functional/funannotate/funannotate/annotate_results/Gryllus_firmus.annotations.txt"
]

species = [
    "G.pennsylvanicus",
    "G.firmus"
]

cols_to_keep = [
    "PFAM",
    "InterPro",
    "EggNog",
    "COG",
    "GO Terms"
]

all_dfs = []

for file_path, sp in zip(files, species):

    df = pd.read_csv(
        file_path,
        sep="\t",
        usecols=cols_to_keep,
        dtype=str,
        na_values=["", " "]
    )

    total_genes = len(df)

    summary = (
        df.notna()
          .sum()
          .reset_index()
    )
    
    summary.columns = ["Annotation_Identifier", "Count"]
    
    summary["Total"] = total_genes
    summary["Percentage"] = (summary["Count"] / total_genes) * 100
    summary["Species"] = sp
    
    # Keep everything
    all_dfs.append(
        summary[
            ["Annotation_Identifier", "Species", "Count", "Total", "Percentage"]
        ]
    )


# -----------------------
# LONG FORMAT (for TSV)
# -----------------------

long_df = pd.concat(all_dfs, ignore_index=True)

long_df.to_csv("annotation_summary_long.tsv", sep="\t", index=False)

# -----------------------
# WIDE FORMAT (for LaTeX)
# -----------------------

wide_df = long_df.pivot(
    index="Annotation_Identifier",
    columns="Species",
    values="Percentage"
).reset_index()

# Format percentages with %
for col in species:
    wide_df[col] = wide_df[col].map(lambda x: f"{x:.0f}\\%")

# -----------------------
# Write LaTex table
# -----------------------

latex_table = f"""
\\begin{{table}}[H]
\\centering
\\begin{{tabular}}{{lcc}}
    \\toprule
    \\textbf{{Annotation Identifier}} & \\textit{{{species[0]}}} & \\textit{{{species[1]}}} \\\\
    \\midrule
"""

for _, row in wide_df.iterrows():
    latex_table += f"    {row['Annotation_Identifier']} & {row[species[0]]} & {row[species[1]]} \\\\\n"

latex_table += f"""    \\bottomrule
\\end{{tabular}}
\\caption{{Genome functional annotation assignment sources for \\textit{{{species[1]}}} and \\textit{{{species[0]}}} annotated gene sets.}}
\\label{{tab:functional_annotation}}
\\end{{table}}
"""

with open("annotation_summary.tex", "w") as f:
    f.write(latex_table)

print("✔ Long TSV saved as annotation_summary_long.tsv")
print("✔ LaTeX table saved as annotation_summary.tex")
