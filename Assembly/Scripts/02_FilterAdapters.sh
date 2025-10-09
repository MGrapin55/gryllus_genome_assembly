#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=HifiAdapterFilter
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks=1                          
#SBATCH --cpus-per-task=8
#SBATCH --mem=16gb
#SBATCH --partition=batch
set -eou pipefail


# ENVIRONMENT SETUP (SWAN HCC)
module purge
module load blast/2.15
module load bamtools/2.5

# ====================================================================================================== #
# Purpose: Uses HifiAdapterFilt.sh to filter adapters from hifi reads 
# ====================================================================================================== #

# Central Location
WKDIR=""

# File prefix (no extension included)
PRE="hifi_reads" 

# Clone the repository 'git clone https://github.com/sheinasim-USDA/HiFiAdapterFilt'
# Path to HiFiAdapterFilt directory
HiFiAdapterFilt=""


##########################################################################################################
##											Commands    												##
##########################################################################################################
# Add HiFiAdapterFilt and DB to PATH
export PATH="$PATH:$HiFiAdapterFilt:$HiFiAdapterFilt/DB"

# Path to raw HiFi BAM file
data_dir="$WKDIR/00_Raw_Data"

# Final output directory on /work
output_dir="$WKDIR/02_Filtered_Adaptors"

################################################################################
##                             Scratch Directory                              ##
################################################################################

scratch_dir="/scratch/$USER/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
work_dir="$scratch_dir/$PRE"  # Where HiFiAdapterFilt will run
mkdir -p "$work_dir"
cd "$work_dir"

echo "Working in scratch directory: $(pwd)"

################################################################################
##                              Copy Input File                              ##
################################################################################

src_file="$data_dir/${PRE}.bam"

if [[ ! -f "$src_file" ]]; then
  echo "ERROR: Input file not found: $src_file" >&2
  exit 1
fi

cp "$src_file" .

echo "Copied $(basename "$src_file") to $(pwd)"

################################################################################
##                          Run HiFiAdapterFilt                              ##
################################################################################

bash "$HiFiAdapterFilt/hifiadapterfilt.sh" \
  -p "$PRE" \
  -t "$SLURM_CPUS_PER_TASK" \
  -o .

################################################################################
##                            Copy Results Back                               ##
################################################################################

mkdir -p "$output_dir"

for suffix in contaminant.blastout blocklist filt.fastq.gz stats; do
  f="${PRE}.${suffix}"
  if [[ -f "$f" ]]; then
    rsync -ah --progress "$f" "$output_dir/"
  else
    echo "WARNING: Expected output file missing: $f" >&2
  fi
done

################################################################################
##                            Reporting & Cleanup                             ##
################################################################################

echo "Done Running $SLURM_JOB_NAME"
echo "Results synced to: $output_dir"
echo "Proceed to SCRIPTS/03_CountKmers"
