#!/bin/bash
#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --job-name=busco
#SBATCH --error=WKDIR/std.err/%x.%J.err
#SBATCH --output=WKDIR/std.out/%x.%J.out
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=8                   
#SBATCH --mem=25gb                           
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load anaconda/25.3

set +eou pipefail
source ~/.bashrc
conda activate "$NRDSTOR/busco"
set -eou pipefail

# Usage: busco_genome.sh <fasta> <lineage> <WKDIR>

# Input arguments
FASTA="$1"
LINEAGE="$2"
WKDIR="$3"

echo "Running BUSCO..."
echo "  FASTA:   $FASTA"
echo "  Lineage: $LINEAGE"
echo "  WKDIR:   $WKDIR"

# Make BUSCO output directory
OUT="$WKDIR/busco_$(basename "$LINEAGE")"
mkdir -p "$OUT"

# Extract FASTA base name without common extensions
NAME=$(basename "$FASTA")
NAME="${NAME%.fa}"
NAME="${NAME%.fasta}"
NAME="${NAME%.fna}"
NAME="${NAME%.fai}"

# Run BUSCO
busco \
    -i "$FASTA" \
    -m genome \
    -l "$LINEAGE" \
    -c "${SLURM_CPUS_PER_TASK}" \
    -o "$NAME" \
    --miniprot \
    --force \
    --out_path "$OUT"