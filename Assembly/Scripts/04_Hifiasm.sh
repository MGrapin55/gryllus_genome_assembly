#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=HiFiasm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=100gb		
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=WKDIR/std.err/%x.%J.err
#SBATCH --output=WKDIR/std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP (HCC SWAN)
module purge
module load hifiasm/0.25

# ====================================================================================================== #
# Purpose: Uses Hifiasm to generate initial phased and hapolype resolved genome assembly 
# ====================================================================================================== #

# Working Directory For Analysis
WKDIR=""

# Path to Filter Reads
FILTERED_READS=""

# Prefix to Name you Inital Assembly 
PREFIX=""


THREADS="$SLURM_CPUS_PER_TASK"

################################################################################
#                            Run HiFiAsm Assembly                              #
################################################################################
# Scratch Directory
SCR="/scratch/$USER_${SLURM_JOB_NAME}"
echo "Scratch directory: $SCR"

OUTDIR=$WKDIR/Hifiasm

# Copy Filtered Reads
[ -e "$FILTERED_READS" ] || echo "$(basename "$FILTERED_READS") does not exist."
cp "$FILTERED_READS" $SCR
echo "Copied $(basename $FILTERED_READS) to $SCR"

echo "Running HiFiAsm..."
hifiasm -t "$THREADS" -o "$PREFIX.asm" "$SCR/$(basename "$FILTERED_READS")"

# Sync Results Back
echo "$(date): Syncing per-run outputs to: $OUTDIR"
rsync -au --exclude "$(basename FILTERED_READS)" $SCR $OUTDIR

echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"
echo 'Proceed to SCRIPTS/05_Assembly_QC'
