#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=blastn
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32                  
#SBATCH --mem=55gb                           
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch

# This script woukd be better if it was converted to a slurm array.

set -euo pipefail
# Might have to play with parametes for the production run 
# adjust for production 32
# adjust for production 110
# ===================== MODULES ===================== #
module purge
module load biodata/1.0
module load blast/2.15
module load anaconda/25.3

set +euo pipefail
source ~/.bashrc
conda activate "$NRDSTOR/biostar"
set -euo pipefail

# ===================== CONFIG ===================== #
WKDIR=""
FASTA=""
BLAST_DB="/work/HCC/BCRF/BLAST"

TMPDIR="/scratch/${USER}_blast"
DBDIR="$TMPDIR/db"
CHUNKDIR="$TMPDIR/chunks"
RESULTSDIR="$TMPDIR/results"
mkdir -p "$TMPDIR/db" "$TMPDIR/chunks" "$TMPDIR/results"

THREADS_PER_JOB=8           # Threads per BLAST job
DB_RAM_GB=10                # Estimated DB memory footprint (initially 25)
PER_THREAD_OVERHEAD_GB=0.15

# ===================== ARGUMENT PARSING ===================== #
NAME=$(basename "$FASTA" | sed 's/\.[^.]*$//')
START_CHUNK=1
while [[ $# -gt 0 ]]; do
  case $1 in
    --start-chunk) START_CHUNK="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

PROGRESS_FILE="$WKDIR/blast_progress_${SLURM_JOB_ID}.txt"
touch "$PROGRESS_FILE"

# ===================== WINDOWMASKER ===================== #
echo "[STEP 1] Running windowmasker..."
# Run windowmasker to build counts/statistics
windowmasker -in "$FASTA" \
             -infmt fasta \
             -mk_counts \
             -sformat obinary \
             -out $WKDIR/${NAME}.counts

# Run windowmasker to mask using counts/statistics
windowmasker -in "$FASTA" \
             -infmt fasta \
             -ustat ${NAME}.counts \
             -dust T \
             -outfmt fasta \
             -out $WKDIR/${NAME}.masked.fa

echo "[DONE WINDOWMASKER]"
# ===================== PREP DATABASE ===================== #
echo "[STEP 2] Copying BLAST DB..."
parallel -j "$SLURM_CPUS_PER_TASK" rsync -ah {} "$TMPDIR/db/" ::: "$BLAST_DB"/nt.*
parallel -j 32 rsync -ah {} "$TMPDIR/db/" ::: "$BLAST_DB"/nt.*

# Copy masked fasta
rsync -a $WKDIR/${NAME}.masked.fa $TMPDIR

echo "[DONE PREPARING DATABASE]"

# ===================== SPLIT FASTA: 1 sequence per file ===================== #ls
echo "[STEP 3] Splitting masked FASTA into individual sequences..."
cd "$TMPDIR/chunks"
awk '/^>/{if(out){close(out)} out=sprintf("seq_%06d.fa", ++i)} {print > out}' $TMPDIR/${NAME}.masked.fa

echo "[DONE SPILITTING FASTA]"

TOTAL_SEQS=$(ls "$TMPDIR/chunks"/*.fa | wc -l)
CPUS=$SLURM_CPUS_PER_TASK
RAM_GB=$(( ${SLURM_MEM_PER_NODE//[!0-9]/} - 5 ))
MEM_PER_JOB_GB=$(awk -v db="$DB_RAM_GB" -v th="$THREADS_PER_JOB" -v oh="$PER_THREAD_OVERHEAD_GB" \
  'BEGIN { printf "%.2f", db + (th * oh) }')
MAX_JOBS=$(awk -v ram="$RAM_GB" -v memjob="$MEM_PER_JOB_GB" 'BEGIN { print int(ram / memjob) }')
[[ $MAX_JOBS -gt $(( CPUS / THREADS_PER_JOB )) ]] && MAX_JOBS=$(( CPUS / THREADS_PER_JOB ))
[[ $MAX_JOBS -lt 1 ]] && MAX_JOBS=1

echo "  "
echo "##################################################################"
echo "[INFO] Total sequences: $TOTAL_SEQS, Max parallel jobs: $MAX_JOBS"
echo "##################################################################"
echo "  "
# ===================== BLAST FUNCTION ===================== #
run_chunk() {
    chunkfile="$1"
    blastn -db "$TMPDIR/db/nt" \
           -query "$chunkfile" \
           -outfmt '6 qseqid staxids bitscore std' \
           -max_target_seqs 10 \
           -max_hsps 1 \
           -evalue 1e-25 \
           -num_threads "$THREADS_PER_JOB" \
           -out "$TMPDIR/results/$(basename "$chunkfile" .fa).blastn.out"
    echo "$(basename "$chunkfile")" >> "$PROGRESS_FILE"
}
# Exporting Variables
export -f run_chunk
export TMPDIR THREADS_PER_JOB PROGRESS_FILE
export DBDIR CHUNKDIR RESULTSDIR WKDIR NAME

# ===================== CLEANUP FUNCTION ===================== #
cleanup() {
    echo "[CLEANUP] EXIT trap triggered."

    # Relax error handling inside trap so partial failures don't abort early
    set +e

    # Check how many chunks finished
    TOTAL_CHUNKS=$(ls "$CHUNKDIR"/seq_*.fa 2>/dev/null | wc -l)
    DONE_CHUNKS=$(sort -u "$PROGRESS_FILE" 2>/dev/null | wc -l)

    if [[ $DONE_CHUNKS -lt $TOTAL_CHUNKS ]]; then
        NEXT_CHUNK=$(( DONE_CHUNKS + 1 ))
        echo "================================================="
        echo " Job ended before finishing all sequences."
        echo " Last completed sequence: $DONE_CHUNKS"
        echo " To restart, run:"
        echo "   sbatch blastn.sh --start-chunk $NEXT_CHUNK"
        echo "================================================="
    fi

    # Combine results if they exist
    if ls "$RESULTSDIR"/*.blastn.out >/dev/null 2>&1; then
        echo "[STEP 5] Combining results..."
        cat "$RESULTSDIR"/*.blastn.out > $WKDIR/${NAME}.ncbi.blastn.out
        echo "BLAST complete. Output saved to: $WKDIR/${NAME}.ncbi.blastn.out"
    else
        echo "[WARNING] No BLAST result files found in $RESULTSDIR"
    fi

    # Copy everything back to WKDIR for inspection/restart
    echo "[CLEANUP] Copying scratch data back to $WKDIR"
    mkdir -p $WKDIR/results_${SLURM_JOB_ID}
    rsync -a "$RESULTSDIR/" $WKDIR/results_${SLURM_JOB_ID}/

    echo "[CLEANUP] Done."
}

# ===================== RUN IN PARALLEL ===================== #
cd "$TMPDIR/chunks"
CHUNK_LIST=$(ls seq_*.fa | sort | tail -n +$START_CHUNK)

echo "[STEP 4] Running BLAST for $(echo $CHUNK_LIST | wc -w) sequences..."
parallel --jobs "$MAX_JOBS" run_chunk ::: $CHUNK_LIST

trap cleanup EXIT

