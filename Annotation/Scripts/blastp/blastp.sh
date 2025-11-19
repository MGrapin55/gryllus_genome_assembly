#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=functionalBlastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch


set -euo pipefail

# ENV (HCC SWAN)
module purge 
module load blast/2.15


# Requirments: (Check for most recent versions)
# NCBI BLAST SUIT

# Michael Grapin @ Moore Lab Research Technician 
# October 8th 2025
####################################################################################################################
##											Parameters and Setup												  ##
####################################################################################################################
# Fill out required fields and then submit via slurm 
# Usage sbatch Blastp.sh 


# Path to where you want to output the breaker3 run (Change if you don't want it in the same directory that you submit the script
OUTDIR=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP

# Braker Amino Acid File
brakerAA=/lustre/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Braker3/Braker3_Run/braker.aa

# Your GFF3 file (from funannotate)
GFF3=/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/funannotate/funannotate/annotate_results/Gryllus_pennsylvanicus.gff3

#SCR=/scratch/$USER/blastp

# Copy files back on exit
# trap 'rsync -auv "/scratch/mgrapin2/blastp" "/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/Gpenn/Gpenn_Production/ANNOTATION/Functional/BlastP"' exit

####################################################################################################################
##												Commands														  ##
####################################################################################################################

# Make directories
mkdir -p $OUTDIR 
#mkdir -p $SCR
cd $OUTDIR

# Copy to scratch
# rsync -au $OUTDIR $SCR
# rsync -a $brakerAA $SCR
# rync -a $GFF3 $SCR
# cd $SCR

# Make BLAST database
# Swissprot Insecta 10,058 proteins
makeblastdb -in swissprot_insecta.fasta -dbtype prot -input_type fasta -out swissprot_insecta

# TrEMBL Insecta 7,071,771 proteins
makeblastdb -in unreviewed_insecta.fasta -dbtype prot -input_type fasta -out unreviewed_insecta

# Swissprot, no taxononomy filtering 573,661 proteins
makeblastdb -in uniprot_sprot.fasta -dbtype prot -input_type fasta -out uniprot_sprot

# Run BLASTp
# Swissprot Insecta
blastp -num_threads $SLURM_CPUS_PER_TASK -query $brakerAA -db swissprot_insecta -outfmt '6 std stitle' -out $OUTDIR/swissprot_insecta_results.blastp -evalue 1e-6 \
-max_hsps 1 -max_target_seqs 1

# TrEMBL Insecta
blastp -num_threads $SLURM_CPUS_PER_TASK -query $brakerAA -db unreviewed_insecta -outfmt '6 std stitle' -out $OUTDIR/unreviewed_insecta_results.blastp -evalue 1e-6 \
-max_hsps 1 -max_target_seqs 1

# Swissprot no taxonomy
blastp -num_threads $SLURM_CPUS_PER_TASK -query $brakerAA -db uniprot_sprot -outfmt '6 std stitle' -out $OUTDIR/uniprot_sprot_results.blastp -evalue 1e-6 \
-max_hsps 1 -max_target_seqs 1


# Pull IDs of transcripts with product "hypothetical protein"
cat $GFF3 | grep "product=hypothetical protein" | cut -f9 | cut -f1 -d ";" | sed s/ID=//g > hypothetical_proteins.txt

echo "[Number of hypothetical proteins in $GFF3]"
echo "$(wc -l hypothetical_proteins.txt)"

echo "######################################################################"
echo "Done Running blastp.sh -> Run SearchBlastResults.sh in $OUTDIR"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"


