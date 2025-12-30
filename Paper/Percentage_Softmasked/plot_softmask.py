import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load your TSV file
df = pd.read_csv("softmask_stats.tsv", sep="\t")

# Make Region a categorical variable and order it properly
region_order = [f"Chromosome {i}" for i in range(1, 16)] + ["Unplaced Scaffolds"]
df["Region"] = pd.Categorical(df["Region"], categories=region_order, ordered=True)

# Rename Chromosome 1 to X Chromosome 
df["Region"] = df["Region"].replace("Chromosome 1", "X Chromosome")

# Rename the Species to sorthand names 
mapping = {
    "Gpenn.chr.final.masked.numt.sorted.renamed": "G.pennsylvanicus",
    "Gfirm.chr.final.masked.numt.sorted.renamed": "G.firmus",
    "GCA_046254815.1_ASM4625481v1_genomic": "G.assimilis",
    "GCA_965638035.1_iqGryBima1.hap1.1_genomic.fna": "G.bimaculatus"
}

df["Species"] = df["Species"].replace(mapping)

# Make a color pallete for the different species 
colors = {
    "G.pennsylvanicus": "#618B4A",
    "G.firmus": "#AFBC88",
    "G.assimilis": "#7AA095",
    "G.bimaculatus": "#49306B"
}

plt.figure(figsize=(14, 6))

sns.barplot(
    data=df,
    x="Region",
    y="Softmask_%",
    hue="Species",
    dodge=True,        # species side-by-side
    linewidth=1,       # edge width of bars (optional)
    edgecolor="black",  # make edges visible (optional)
    palette=colors
)


plt.xticks(rotation=45, ha="right")
plt.xlabel("")
plt.ylabel("Softmask (%)")
plt.title("Softmask Percentage by Region Across Species")
plt.tight_layout()

# Save as PNG
plt.savefig("softmask_boxplot.png", dpi=300)  # high-resolution PNG

# Optionally show the plot
plt.show()
