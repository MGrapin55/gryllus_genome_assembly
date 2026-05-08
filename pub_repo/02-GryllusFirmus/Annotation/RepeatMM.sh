#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=Gfirm_RepeatMasker
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load anaconda
conda activate $NRDSTOR/repeat_annnotation

# Michael Grapin @ Moore Lab UNL 
# ====================================================================================================== #
# Purpose: Annotate Repeat Regions of a Genome Assembly
# ====================================================================================================== #

# [Altered Script to just run the RepeatMasker Step]
WKDIR=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/RepeatMM/Masked

FASTA=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/scaffolding/Ragtag_assimilis_parameters/final/Gfirm.chr.fa

LIB=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gfirm/Gfirm_Production/ANNOTATION/RepeatMM/RepeatLibrary/Gfirm_families.prefix.fa.known.unknown.FINAL

###########################################################################################################
##                                           Commamds                                                    ##
###########################################################################################################
cd $WKDIR
#run repeatmasker
RepeatMasker -pa $SLURM_CPUS_PER_TASK -gff -s -a -inv -no_is -norna -xsmall -nolow -div 40 -lib $LIB -cutoff 225 $FASTA


echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"
