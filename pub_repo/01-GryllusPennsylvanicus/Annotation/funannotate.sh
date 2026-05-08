#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=Funannotate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
#SBATCH --qos=short

# ENV (HCC SWAN)
module purge 
module load anaconda/25.3
conda activate $NRDSTOR/funannotate

# Moved error catching because unset variable in conda env script
set -euo pipefail

# Requirments: (Check for most recent versions)
# Funannotate (https://anaconda.org/bioconda/funannotate)
# Busco (https://anaconda.org/bioconda/busco)
# AGAT (https://anaconda.org/bioconda/agat)
# GFFTK (https://pypi.org/project/gfftk/)

# Michael Grapin @ Moore Lab Research Technician 
# October 8th 2025
####################################################################################################################
##											Parameters and Setup												  ##
####################################################################################################################
# Fill out required fields and then submit via slurm 
# Usage sbatch funannotate.sh 


# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/funannotate 

# Braker gtf file 
GTF=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Braker3_Run/braker.gtf

# Genome Assembly fasta file 
FASTA=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Gpenn.chr.final.masked.numt.sorted.renamed.fasta

# Eggnog Mapper Annotation Results 
EGGNOG=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/eggNog/out.emapper.annotations

# InterproScan Annotation Results
IPR=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/InterPro/Gpenn.braker.aa.NoStop.IPR.xml

# Busco Lineage (Path or db or just lineage name (ex. insecta)) ## TBD
BUSCO=insecta


# Your Organism's species 
SPECIES=Gpenn



####################################################################################################################
##												Commands														  ##
####################################################################################################################

# Make directories
mkdir -p $OUTDIR $OUTDIR/funannotate
cd $OUTDIR

#echo "Converting Braker GFF..."
#GFF3=$(basename $GTF .gtf)
# Convert braker output using agat
#agat_convert_sp_gxf2gxf.pl -g $GTF -o $OUTDIR/${SPECIES}_${GFF3}.gff3

#echo "Sanitizing GFF3..."
# Sanitize agat's output
#gfftk sanitize -f $FASTA -g $OUTDIR/${SPECIES}_${GFF3}.gff3 -o $OUTDIR/${SPECIES}_braker.sanitized.gff3

#echo "Check to make sure your gff was converted properly"

echo " Running Funannotate..."
# Run funannotate
funannotate annotate --gff $OUTDIR/${SPECIES}_braker.sanitized.gff3 --fasta $FASTA \
        --eggnog $EGGNOG \
        --iprscan $IPR \
        --busco_db $BUSCO \
        --out $OUTDIR/funannotate \
        --species "Gryllus pennsylvanicus" \
        --cpus $SLURM_CPUS_PER_TASK --force
        
echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"



