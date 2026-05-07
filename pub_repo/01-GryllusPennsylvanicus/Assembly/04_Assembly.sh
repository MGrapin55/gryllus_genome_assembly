#!/bin/bash
#SBATCH --time=167:55:00
#SBATCH --job-name=Hifiasm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --error=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/std.err/%x.%J.err
#SBATCH --output=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load hifiasm/0.25

# ====================================================================================================== #
# Purpose: Uses Hifiasm to generate a genome assembly 
# ====================================================================================================== #

# Central Location
WKDIR="/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn"

FILTERED_READS="hifi_reads.filt.fastq.gz"

RUN_ID="Gpenn_run2"  # <-- CHANGE this for each run (e.g., run2, run3, etc.) (Ex. Gpenn_run1)

################################################################################
#                              Directory Setup                                 #
################################################################################

DATA="$WKDIR/02_Filtered_Adaptors/$FILTERED_READS"
COMMON_DIR="$WKDIR/04_Assembly/common"
RUNS_DIR="$WKDIR/04_Assembly/runs"

PREFIX="${RUN_ID}"
THREADS="$SLURM_CPUS_PER_TASK"

################################################################################
#                          Setup Scratch Directory                             #
################################################################################

scratch_dir="/scratch/$USER/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
run_dir="$scratch_dir/$RUN_ID"
mkdir -p "$run_dir"
cd "$run_dir"

echo "Working in scratch directory: $(pwd)"

################################################################################
#                              Copy Input File                                 #
################################################################################

if [[ ! -f "$DATA" ]]; then
  echo "ERROR: Input FASTQ file not found: $DATA" >&2
  exit 1
fi

cp "$DATA" .
FASTQ="$(basename "$DATA")"

################################################################################
#                        Optionally Copy .bin for Reuse                        #
################################################################################
# removed the prefix
mkdir -p "$COMMON_DIR"
if compgen -G "$COMMON_DIR/${PREFIX%%_*}.asm.*.bin" > /dev/null; then
  echo "Reusing existing .bin files from common directory..."
  cp "$COMMON_DIR/${PREFIX%%_*}.asm."*.bin .
else
  echo "No .bin files found for reuse in $COMMON_DIR."
fi

################################################################################
#                            Run HiFiAsm Assembly                              #
################################################################################

echo "Running HiFiAsm..."
hifiasm -t "$THREADS" -o "$PREFIX.asm" "$FASTQ"

################################################################################
#                            Sync Results Back                                 #
################################################################################

mkdir -p "$RUNS_DIR/$RUN_ID"

echo "$(date): Syncing per-run outputs to: $RUNS_DIR/$RUN_ID"
rsync -au "$run_dir/" "$RUNS_DIR/$RUN_ID/"

# Also update .bin files in the shared common/ directory
echo "$(date): Updating common .bin files in: $COMMON_DIR"
rsync -au --include="*.bin" --exclude="*" "$run_dir/" "$COMMON_DIR/"

################################################################################
#                         Reporting and Optional Cleanup                       #
################################################################################

echo "Done Running $SLURM_JOB_NAME"
echo "Results synced to: $RUNS_DIR/$RUN_ID"
echo "Common .bin files synced to: $COMMON_DIR"

# Optional: cleanup scratch (if desired)
# rm -rf "$scratch_dir"

echo "######################################################################"
echo "=== Memory Report ==="
command -v mem_report &>/dev/null && mem_report || echo "mem_report not found"

echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"
echo 'Proceed to SCRIPTS/05_Assembly_QC'
