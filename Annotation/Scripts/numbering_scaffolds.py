# code from: https://github.com/kataokaklab/Gryllus_assimilis_genome/blob/main/structural_annotation/numbering_scaffolds.py

# Renames sorted genome assemble (fasta) to Super-Scaffold_n...

from Bio import SeqIO
import pandas as pd

### Setting ###
genome = 'genome.sorted.fasta'
output_dir = 'OUTDIR PATH'

### Command ###
fasta_in = open(genome)
counter = 1
for record in SeqIO.parse(fasta_in, 'fasta'):
    seq = record.seq
    with open(f'{output_dir}/Gpenn.chr.final.masked.numt.sorted.renamed.fasta', mode='a') as f:
        f.write(f'>Super-Scaffold_{counter}\n')
        f.write(f'{seq}\n')
    counter += 1
