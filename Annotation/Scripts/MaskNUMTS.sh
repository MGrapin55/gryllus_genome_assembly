#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=MaskNUMTS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb		
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch


set -euo pipefail

# ENV (HCC SWAN)
module purge 
module load blast/2.15
module load bedtools/2.31


# Michael Grapin @ Moore Lab Research Technician 
# September 23th 2025
####################################################################################################################
##											Parameters and Setup												  ##
####################################################################################################################
# Fill out required fields and then submit via slurm 
# Usage sbatch MaskNUMTS.sh 
# Following approach used in Liu et al. (2024, Mol Phyl Evol) (https://doi.org/10.1016/j.ympev.2024.108221)



# Path to your working dirctory
#WKDIR=/scratch/${USER}

# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=

# Your Organism's species 
SPECIES=

# Path to genome assembly
GENOME=

# Path to Mitogenome assembly
MITOGENOME=


####################################################################################################################
##												Commands														  ##
####################################################################################################################
# Make std.err and std.out directory for your breaker run
mkdir -p $OUTDIR

# Move to OUTDIR
cd $OUTDIR 

# Make BLAST database
makeblastdb -in "$GENOME" -dbtype nucl -out "${SPECIES}_clean"

# Run BLASTn
blastn \
 -query $MITOGENOME \
 -db ${SPECIES}_clean \
 -outfmt 6 \
 -perc_identity 80 -evalue 1e-6 \
 -out ${SPECIES}_numt_blast.results \
 -num_threads $SLURM_CPUS_PER_TASK


# Copy Blast Files back to the $OUTDIR
#rsync -au $WKDIR/ $OUTDIR/
#echo "Done Copying Files..." 
#ls -lh $OUTDIR

# Convert Blast results into .bed format
awk '$4 >= 100' ${SPECIES}_numt_blast.results | awk '{start = ($9 < $10 ? $9 : $10); end = ($9 > $10 ? $9 : $10); strand = ($9 < $10 ? "+" : "-"); print $2, start-1, end, strand;}' OFS='\t' > ${SPECIES}_numts.bed

bedtools maskfasta -fi $GENOME -bed ${SPECIES}_numts.bed -fo $GENOME.numt -soft

mv $GENOME.numt $OUTDIR

echo "[Done Masking NUMTs]" 
echo " "
echo " "
echo "[RESOURCE REPORT]"
echo "Run: seff $SLURM_JODID"