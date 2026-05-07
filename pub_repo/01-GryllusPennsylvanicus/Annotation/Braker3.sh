#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=Braker3_run1
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=128gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch


set -euo pipefail

# ENV (HCC SWAN)
module purge 
module load apptainer
export BRAKER_SIF="/home/moorelab/mgrapin2/Apptainer_Images/braker3.sif"
export AUGUSTUS_CONFIG_PATH="/home/moorelab/mgrapin2/Apptainer_Images/Augustus_config"
# Only necessary for older versions of Braker Pipeline 
# export GM_KEY="PATH/TO/GM_KEY"

# Otherwise you need to download the singularity image (.sif) (https://hub.docker.com/r/teambraker/braker3)
# singularity build braker3.sif docker://teambraker/braker3:latest
# singularity exec braker3.sif braker.pl
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test1.sh .
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test2.sh .
# singularity exec -B $PWD:$PWD braker3.sif cp /opt/BRAKER/example/singularity-tests/test3.sh .
# export BRAKER_SIF=/your/path/to/braker3.sif
# Copy the AUGUSTUS_CONFIG
# singularity exec ${BRAKER_SIF} cp -r /opt/Augustus/config ./Augustus_config

# Command in the test#.sh needs a --AUGUSTUS_CONFIG_PATH=<PATH>

# Test the pipeline you wish to run (BRAKER3)
# bash test1.sh # tests BRAKER1
# bash test2.sh # tests BRAKER2
# bash test3.sh # tests BRAKER3


# Michael Grapin @ Moore Lab Research Technician 
# September 18th 2025
####################################################################################################################
##											Parameters and Setup												  ##
####################################################################################################################
# Fill out required fields and then submit via slurm 
# Usage sbatch Braker3.sh 

# Examine /scratch space on runing jon on HCC SWAN 
# srun --jobid=<yourjobid> --pty $SHELL 
# cd /scratch/$USER | ls -lh

# Path to your working dirctory where you want your braker run ( Running on scratch spaces to help speed on in/out operations) 
WKDIR=/scratch/${USER}/Gpenn

# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Braker3_Run

# Set up a trap to copy results back even if job crashes or is canceled
trap 'rsync -auv "/scratch/${USER}/Gpenn/run/" "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Braker3_Run"' exit

# Your Organism's species 
SPECIES=Gpenn

# Path to genome assembly
GENOME=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Gpenn.chr.final.masked.numt.sorted.renamed.fasta

# Path to Bam File(s)
BAM=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/RNAseq/BAM

# Path to your ORTHO_DB protein seqiences (https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)
PROT=/mnt/nrdstor/moorelab/shared/ORTHO_DB/Arthropoda.fa.gz


####################################################################################################################
##												Commands														  ##
####################################################################################################################

# Make directories
mkdir -p $WKDIR $WKDIR/run $OUTDIR 

# Copy inputs to $WKDIR (SCRATCH SPACE)
rsync -au $GENOME $WKDIR
find "$BAM" -name "*.bam" -print0 | xargs -0 -I{} rsync -au {} "$WKDIR/"
rsync -au $PROT $WKDIR

# Move to WKDIR
cd $WKDIR 

# Unzip the $PROT
gunzip $(basename $PROT)

# Get list of BAM files
BAM_LIST=$(find $WKDIR -name "*.bam" | paste -sd, -)

# Run the Breaker Pipline 
singularity exec -B ${WKDIR}:${WKDIR} ${BRAKER_SIF} braker.pl \
--species=$SPECIES \
--genome=$(basename $GENOME) \
--prot_seq=$(basename $PROT .gz) \
--bam $BAM_LIST \
--workingdir=${WKDIR}/run \
--threads $SLURM_NTASKS \
--gff3 \
--AUGUSTUS_CONFIG_PATH=$AUGUSTUS_CONFIG_PATH

# Copy Braker Files back to the $OUTDIR
rsync -au $WKDIR/run/ $OUTDIR/
echo "Done Copying Files..." 
ls -lh $OUTDIR

echo "[Done Braker3 Pipeline]" 
echo " "
echo " "
echo "[RESOURCE REPORT]"
echo "Run: seff $SLURM_JODID"
