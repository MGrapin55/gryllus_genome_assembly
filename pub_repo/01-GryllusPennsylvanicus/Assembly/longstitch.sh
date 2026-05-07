#!/bin/bash
#SBATCH --time=12:00:00          
#SBATCH --job-name=longstitch
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --mail-user=
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
conda activate "$NRDSTOR/longstitch"
set -eou pipefail

# ====================================================================================================== #
# Purpose: Run longstitch pipeline to scaffold contig assembly
# ====================================================================================================== #

# Central Location
WKDIR="/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/scaffolding"

# Genome size (bp) # chosen from previous scaffolding 
size=2185179992

# Draft Assembly (fasta format)
draft="/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/blob_fasta/mitogenome_removed/Gpenn.p.clean.mito.fa"

# Raw Reads
reads="/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/02_Filtered_Adaptors/hifi_reads.filt.fastq.gz"

# Prefix for scaffolds
Prefix="Gpenn"

# Version
V=default
##########################################################################################################
##											Commands    												##
##########################################################################################################
mkdir -p $WKDIR/longstitch.V${V}
cd $WKDIR/longstitch.V${V}

# Check if file exists
[ -e "$draft" ] || { echo "$draft doesn't exist."; exit 1; }
[ -e "$reads" ] || { echo "$reads doesn't exist."; exit 1; }

# get file name
d=$(basename "$draft")
r=$(basename "$reads")
# create soft link
ln -sf "$draft" "$d"
ln -sf "$reads" "$r"

echo "# Directory Contents"
echo "# ----------------------------------------------------------------------------------- #"
ls -l
echo "# ----------------------------------------------------------------------------------- #"

# longstitch pipeline command 
echo "Starting Longstitch Run..."
longstitch ntLink-arks draft="${d%.fa}" reads="${r%.fastq.gz}" G="$size" longmap=hifi out_prefix="$Prefix" t=$SLURM_CPUS_PER_TASK rounds=5 gap_fill=True

echo "Done with longstitch!"
