#!/bin/bash
#SBATCH --time=24:00:00          # Run time in hh:mm:ss
#SBATCH --job-name=CountKmers
#SBATCH --error=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/std.err/%x.%J.err
#SBATCH --output=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/std.out/%x.%J.out
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1                          
#SBATCH --cpus-per-task=16
#SBATCH --mem=50gb                           # trial with 50gb
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load meryl/1.4

# ====================================================================================================== #
# Purpose: Uses Meryl to count Kmers from Hifi Reads 
# ====================================================================================================== #

# Central Location
WKDIR="/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn"


# Argument: k-mer size (Fill in at command line) (can use $1 at command line)
KMER=21

# Fastq Name (Update later to handle for then one)
FASTQ="hifi_reads.filt.fastq.gz"

##########################################################################################################
##											Commands    												##
##########################################################################################################
POST_SCRIPT="$WKDIR/SCRIPTS/POST_MERYL.sh"

# Input FASTQ (gzipped OK)
FASTA1="$WKDIR/02_Filtered_Adaptors/$FASTQ"

# Permanent output directory
OUTDIR="$WKDIR/03_Count_Kmers"         
mkdir -p "$OUTDIR"

# Scratch directory
SCRATCH_DIR="/scratch/$USER/$SLURM_JOB_ID/${KMER}mer_Files"
mkdir -p "$SCRATCH_DIR"


# Output DB prefix
DB_NAME="${KMER}mer_db.meryl"

##########################################################################################################
##											Run Meryl													##
##########################################################################################################

echo "$(date): Counting ${KMER}-mers with Meryl..."
# Might Change the <memory=n> and not include it 
meryl count \
    k=$KMER \
    output $SCRATCH_DIR/$DB_NAME \
    $FASTA1


echo "$(date): Done counting ${KMER}-mers."

##########################################################################################################
##										Rsync to Permanent Location										##
##########################################################################################################

echo "$(date): Syncing results to $OUTDIR using rsync..."
rsync -au "$SCRATCH_DIR" "$OUTDIR"

echo "$(date): Done with ${KMER}-mers."
echo "Output synced to: $OUTDIR"

##########################################################################################################
##											Sumbit Meryl Reports										##
##########################################################################################################

# Submit report job with correct paths and dependency
PERM_DIR="$OUTDIR/${KMER}mer_Files"

POST_JOB_ID=$(sbatch --dependency=afterok:$SLURM_JOB_ID --parsable \
    "$POST_SCRIPT" $KMER "$PERM_DIR")

echo "Submitted post-processing job: $POST_JOB_ID"

##########################################################################################################
##											Job Wrap-Up													##
##########################################################################################################
# Where to write the log
LOGFILE="$OUTDIR/${KMER}mer_Files/run_log.txt"

{
    echo "==========================================="
    echo "Meryl Count + Report Job Log"
    echo "Date: $(date)"
    echo "K-mer Size: $KMER"
    echo "Count Job ID: $SLURM_JOB_ID"
    echo "Check with: seff $SLURM_JOB_ID"
    echo "Post-Processing Job ID: $POST_JOB_ID"
    echo "Check with: seff $POST_JOB_ID"
    echo "Output Directory: $OUTDIR/${KMER}mer_Files"
    echo "==========================================="
} >> "$LOGFILE"


echo "######################################################################"
echo "Kmer Size: $KMER"
echo "=== Memory Report ==="
mem_report || echo "mem_report not found"

echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"
echo "Job Finished at: $(date)"
echo "Proceed to SCRIPTS/04_Assembly.sh"