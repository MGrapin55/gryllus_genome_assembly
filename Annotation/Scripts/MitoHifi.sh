#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=MitoHifi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64gb		
#SBATCH --mail-user=mgrapin2@nebraska.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --error=WKDIR/std.err/%x.%J.err
#SBATCH --output=WKDIR/std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load apptainer

# Install Command (Fill in as necessary) 
# apptainer pull mitohifi_<Version>.sig docker://ghcr.io/marcelauliano/mitohifi:master
# Github: https://github.com/marcelauliano/MitoHiFi/tree/master


# ====================================================================================================== #
# Purpose: Extract and Annotatoe the mitochrondiral genomem with the MitoHifi Pipeline 
# ====================================================================================================== #

# Path to MitoHifi Docker Image (.sif)
MITO="/home/moorelab/mgrapin2/Apptainer_Images/mitohifi_3.2.1.sif"

# Output Directory
OUTDIR=

# Version of MitoGenome 
V=

# Assembled fasta contigs/scaffolds to be searched to find mitogenome
FASTA=

# Closely releated ancestor mitogenome in Fasta Format (.fasta, .fa, .fna, etc...)
ANCESTOR_FA=

# Closely releated ancestor mitogenome in GeneBank Format (.gbk)
ANCESTOR_GBK=

# Annotaion format of Genetic Code
CODE=5

#   -o <GENETIC CODE>     -o: Organism genetic code following NCBI table (for
#                           mitogenome annotation): 
#                             1. The Standard Code 
#                             2. The Vertebrate Mitochondrial Code 
#                             3. The Yeast Mitochondrial Code 
#                             4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the
#                                Mycoplasma/Spiroplasma Code 
#                             5. The Invertebrate Mitochondrial Code 
#                             6. The Ciliate, Dasycladacean and Hexamita Nuclear Code 
#                             9. The Echinoderm and Flatworm Mitochondrial Code 
#                            10. The Euplotid Nuclear Code 
#                            11. The Bacterial, Archaeal and Plant Plastid Code 
#                            12. The Alternative Yeast Nuclear Code 
#                            13. The Ascidian Mitochondrial Code 
#                            14. The Alternative Flatworm Mitochondrial Code 
#                            16. Chlorophycean Mitochondrial Code 
#                            21. Trematode Mitochondrial Code 
#                            22. Scenedesmus obliquus Mitochondrial Code 
#                            23. Thraustochytrium Mitochondrial Code 
#                            24. Pterobranchia Mitochondrial Code 
#                            25. Candidate Division SR1 and Gracilibacteria Code

############################################################################################################
##                                              Commands                                                  ##
############################################################################################################
mkdir -p $OUTDIR
cd $OUTDIR

echo "Starting the MitoHifi Pipeline..."
apptainer exec $MITO mitohifi.py -c $FASTA -f $ANCESTOR_FA -g $ANCESTOR_GBK -t $SLURM_CPUS_PER_TASK -o $CODE

echo "Done MitoHifi Pipeline..."
echo "Check $OUTDIR for results..."

# usage: MitoHiFi [-h] (-r <reads>.fasta | -c <contigs>.fasta) -f
#                 <relatedMito>.fasta -g <relatedMito>.gbk -t <THREADS> [-d]
#                 [-a {animal,plant,fungi}] [-p <PERC>] [-m <BLOOM FILTER>]
#                 [--max-read-len MAX_READ_LEN] [--mitos]
#                 [--circular-size CIRCULAR_SIZE]
#                 [--circular-offset CIRCULAR_OFFSET] [-winSize WINSIZE]
#                 [-covMap COVMAP] [-v] [-o <GENETIC CODE>]
#
# required arguments:
#   -r <reads>.fasta      -r: Pacbio Hifi Reads from your species
#   -c <contigs>.fasta    -c: Assembled fasta contigs/scaffolds to be searched
#                           to find mitogenome
#   -f <relatedMito>.fasta
#                           -f: Close-related Mitogenome is fasta format
#   -g <relatedMito>.gbk  -k: Close-related species Mitogenome in genebank
#                           format
#   -t <THREADS>          -t: Number of threads for (i) hifiasm and (ii) the
#                           blast search
#
# optional arguments:
#   -d                    -d: debug mode to output additional info on log
#   -a {animal,plant,fungi}
#                           -a: Choose between animal (default) or plant
#   -p <PERC>             -p: Percentage of query in the blast match with close-
#                           related mito
#   -m <BLOOM FILTER>     -m: Number of bits for HiFiasm bloom filter [it maps
#                           to -f in HiFiasm] (default = 0)
#   --max-read-len MAX_READ_LEN
#                           Maximum lenght of read relative to related mito
#                           (default = 1.0x related mito length)
#   --mitos               Use MITOS2 for annotation (opposed to default
#                           MitoFinder
#   --circular-size CIRCULAR_SIZE
#                           Size to consider when checking for circularization
#   --circular-offset CIRCULAR_OFFSET
#                           Offset from start and finish to consider when looking
#                           for circularization
#   -winSize WINSIZE      Size of windows to calculate coverage over the
#                           final_mitogenom
#   -covMap COVMAP        Minimum mapping quality to filter reads when building
#                           final coverage plot
#   -v, --version         show program's version number and exit
#   -o <GENETIC CODE>     -o: Organism genetic code following NCBI table (for
#                           mitogenome annotation): 
#                             1. The Standard Code 
#                             2. The Vertebrate Mitochondrial Code 
#                             3. The Yeast Mitochondrial Code 
#                             4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the
#                                Mycoplasma/Spiroplasma Code 
#                             5. The Invertebrate Mitochondrial Code 
#                             6. The Ciliate, Dasycladacean and Hexamita Nuclear Code 
#                             9. The Echinoderm and Flatworm Mitochondrial Code 
#                            10. The Euplotid Nuclear Code 
#                            11. The Bacterial, Archaeal and Plant Plastid Code 
#                            12. The Alternative Yeast Nuclear Code 
#                            13. The Ascidian Mitochondrial Code 
#                            14. The Alternative Flatworm Mitochondrial Code 
#                            16. Chlorophycean Mitochondrial Code 
#                            21. Trematode Mitochondrial Code 
#                            22. Scenedesmus obliquus Mitochondrial Code 
#                            23. Thraustochytrium Mitochondrial Code 
#                            24. Pterobranchia Mitochondrial Code 
#                            25. Candidate Division SR1 and Gracilibacteria Code
echo "######################################################################"
echo "=== Memory Report ==="
command -v mem_report &>/dev/null && mem_report || echo "mem_report not found"

echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"