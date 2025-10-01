#!/bin/bash
#SBATCH --time=168:00:00          # Run time in hh:mm:ss
#SBATCH --job-name=quast
#SBATCH --error=WKDIR/std.err/%x.%J.err
#SBATCH --output=WKDIR/std.out/%x.%J.out
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1                         
#SBATCH --cpus-per-task=8                   
#SBATCH --mem=50gb                           
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load quast/5.0

# Usage: quastL.sh <fasta> <ref> <label> <WKDIR>

FASTA="$1"
REF="$2"
LABEL="$3"
WKDIR="$4"

echo "Running QUAST..."
echo "  FASTA:  $FASTA"
echo "  Ref:    $REF"
echo "  Label:  ${LABEL:-<none>}"
echo "  WKDIR:  $WKDIR"

# Make Quast output directory
OUT="$WKDIR/quast"
mkdir -p "$OUT"

# Example command (adjust as needed)
if [[ -n "$REF" ]]; then
    quast.py "$FASTA" -r "$REF" ${LABEL:+--label "$LABEL"} -o "$OUT" -t "$SLURM_CPUS_PER_TASK"
else
    quast.py "$FASTA" ${LABEL:+--label "$LABEL"} -o "$WKDIR/quast_out" -t "$SLURM_CPUS_PER_TASK"
fi