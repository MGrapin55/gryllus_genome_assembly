#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --job-name=HPC_Blastp_redo
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
HPC1=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl

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
# [First!] Have all commands run on the hpc_FASTA_GridRunner.pl

# Fill out required fields and then submit via slurm
# Usage sbatch redo_jobs.sh 

# HPC GRID RUNNER Config File 
CONFIG=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC/config.conf

# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC

# Unfinished HPC Grid Runner Commands 
CMD=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP/HPC/braker_unreviewed_insecta_blastp_search.cmds.failures

####################################################################################################################
##												Commands														  ##
####################################################################################################################

# Run BLASTp
$HPC1 \
-c $CMD \
-G $CONFIG \

# Find and Concatentate Output
echo "# After All Jobs Complete Run:" 
#echo "find $OUTDIR/${query_no_ext}_${db_no_ext}_blastp_search -name "*.fa.OUT" -exec cat {} \; > $OUTDIR/${query_no_ext}_${db_no_ext}_all.blast.out"


echo "######################################################################"
echo "Done Running blastp.sh -> Run SearchBlastResults.sh in $OUTDIR"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"


