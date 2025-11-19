#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=HPC_Blastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch


set -euo pipefail

# ENV (HCC SWAN)
module purge 
module load blast/2.15
# Path to HPC Grid Runner Fasta Script 
HPC=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC/HpcGridRunner-1.0.2/BioIfx/hpc_FASTA_GridRunner.pl

# [IMPORTANT] 
# CHANGE HpcGridRunner-1.0.2/BioIfx/../PerlLib/HPC/FarmIt.pm line 91 to my $host = hostname() || $ENV{HOSTNAME};

# Requirments: (Check for most recent versions)
# NCBI BLAST SUIT
# HPC GRID RUNNER

# Michael Grapin @ Moore Lab Research Technician 
# November 7th 2025
####################################################################################################################
##											Parameters and Setup												  ##
####################################################################################################################
# [First!] Make a blast db

# Fill out required fields and then submit via slurm
# Usage sbatch blastp_grid.sh 

# HPC GRID RUNNER Config File 
CONFIG=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC/config.conf

# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC

# Query Sequences (Braker Amino Acid File)
QUERY=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Braker3_Run/braker.aa

# DB
DB=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/unreviewed_insecta


####################################################################################################################
##												Commands														  ##
####################################################################################################################
query_no_ext="$(basename "${QUERY%.*}")"
db_no_ext="$(basename "${DB%.*}")"

# Run BLASTp
$HPC \
--cmd_template "blastp -num_threads 4 -query __QUERY_FILE__ -db $DB -outfmt '6 std stitle' -evalue 1e-6 -max_hsps 1 -max_target_seqs 1" \
--query_fasta $QUERY \
-G $CONFIG \
-N 10 \
-O $OUTDIR/${query_no_ext}_${db_no_ext}_blastp_search

# Find and Concatentate Output
echo "# After All Jobs Complete Run:" 
echo "find $OUTDIR/${query_no_ext}_${db_no_ext}_blastp_search -name "*.fa.OUT" -exec cat {} \; > $OUTDIR/${query_no_ext}_${db_no_ext}_all.blast.out"


echo "######################################################################"
echo "Done Running blastp.sh -> Run SearchBlastResults.sh in $OUTDIR"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"


