#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=uniprot_download
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch

# ENV (HCC SWAN)
module purge 
module load anaconda/25.3
conda activate $NRDSTOR/prepocess_repeats

echo "Downloading unreviewed fasta..."
python query_uniprot.py

#echo "Downloading swissprot all taxa fasta..."
#curl -o uniprot_sprot.fasta.gz https://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/complete/uniprot_sprot.fasta.gz        


echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"



