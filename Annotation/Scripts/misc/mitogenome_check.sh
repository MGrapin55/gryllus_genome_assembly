#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=mitogenome_consensus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=30gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --error=./%x.%J.err
#SBATCH --output=./%x.%J.out
#SBATCH --partition=batch

set -euo pipefail

module purge
module load minimap2/2.26
module load samtools/1.20
module load mafft/7.526

# Inputs
READS="reads.fastq"
MITO="mito.fasta"
PREFIX="sample1"

# Validate
[[ -f "$READS" ]] || { echo "Missing READS"; exit 1; }
[[ -f "$MITO" ]] || { echo "Missing MITO"; exit 1; }

echo "Running minimap2..."
minimap2 -ax map-hifi -t "$SLURM_CPUS_PER_TASK" "$MITO" "$READS" | \
samtools sort -o "$PREFIX.bam"

samtools index -c "$PREFIX.bam"

echo "Generating Consensus..."
samtools consensus -f fasta "$PREFIX.bam" -o "$PREFIX.consensus.fa" -X hifi

echo "Concatenating sequences..."
cat "$MITO" "$PREFIX.consensus.fa" > "$PREFIX.concat.fa"

echo "Making Alignment..."
mafft --globalpair --maxiterate 1000 --clustalout "$PREFIX.concat.fa" > "$PREFIX.alignment_result.aln"
mafft --globalpair --maxiterate 1000 "$PREFIX.concat.fa" > "$PREFIX.alignment_result.fa"
