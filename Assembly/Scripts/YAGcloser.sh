#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=YAGcloser
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb		
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load minimap2/2.26
module load samtools/1.20
module load anaconda/25.3

set +eou pipefail
source ~/.bashrc
conda activate "$NRDSTOR/yagcloser"
set -eou pipefail

# ====================================================================================================== #
# Purpose: Close gaps with YAGcloser 
# ====================================================================================================== #
# Genome Assembly in Fasta format 
FASTA=

# Raw Reads used to generated the genome assembly
READS=

# Specifify the directory your wish for output
OUTDIR=

# Detgap excutable (https://github.com/dfguan/asset)
DETGAP=

# YagCloser Repo (https://github.com/merlyescalona/yagcloser)
YAG=
###########################################################################################################
##                                           Commamds                                                    ##
###########################################################################################################
mkdir -p $OUTDIR
cd $OUTDIR
NAME=$(basename "$FASTA" | sed 's/\.[^.]*$//')

# Map the long reads against the genome assembly using minimap2 (Li, 2018).
# Convert the alignment into its binary form (BAM file), sort it by coordinates and index it with samtools (Li et al., 2009) [Version tested: 1.7 (using htslib 1.8)]
# Filter a subset of primary alignments with MAPQ > 20 to identify the closable gaps.
# 1. Run minimap2 and save SAM
minimap2 -t 8 -ax map-hifi $FASTA $READS > aln.sam

# 2. Convert to BAM, filter, and sort
samtools view -b -q 20 aln.sam | samtools sort -o aln.s.bam
echo "Done Alignment"

# Index BAM file
samtools index -@ 8 aln.s.bam
echo "Done Indexing"

# Get a description of the gaps present in the reference file, you can use the detgaps tool from asset
$DETGAP $FASTA > gaps.bed
echo " Generated .bed file"

# Run to identify potential gaps/edits that will be done to your reference
python $YAG/yagcloser.py -g $FASTA \
    -a aln.s.bam \
    -b gaps.bed \
    -o $OUTDIR \
    -f 20 -mins 5 -s $NAME
echo " Identified Gaps"
# [REQUIRED]
# -g FASTA FILE PATH, --genome FASTA FILE PATH: Filepath of the reference genome file to be used. Accepts compressed files (GZIP) (default: None)
# -b BED FILE PATH, --bed BED FILE PATH: Filepath of the bed file describing the gaps of the genome. Accepts compressed files (GZIP) (default:None)
# -a BAM FILE PATH, --aln BAM FILE PATH: Filepath of the alignment of reads to the reference genome in BAM format. This file needs to be indexed before running. (default: None)
# -o FOLDER_PATH, --output FOLDER_PATH: Output path folder. (default: None)
# -s STR, --samplename STR: Short sample name that will be used for naming OUTPUT files. (default: None)
    
# [OPTIONAL]
# -mins INT, --min-support INT: Minimum number of reads needed spanning a gap to be considered for closing or filling, (default: 5).
# -f INT, --flanksize INT: Flank size to be used to select the reads that are in the surroundings of the gap and determine whether there are reads that span the gap or not, (default: 20).
# -mapq <MAPQ_threshold>, --mapping-guality-threshold <MAPQ_threshold>: MAPQ value used to filter alignments for posterior processing. Discarding alingments where: alignment_mapq < MAPQ_threshold, (default: 20).
# -mcc INT, --min-coverage-consensus INT: Require that more than INT sequences to be part of an alignment to put in the consensus. Modify if -mins/--min-support < default, (default: 2).
# -prt FLOAT, --percent-reads-threshold FLOAT: Require that more than INT sequences to remain after the length agreement to be considered for consensus. (default: 0.5)
# -eft FLOAT, --empty-flanks-threshold: FLOAT Percentage of empty flanks required to skip an ambiguous decision on a gap. (default: 0.2)
# -l <log_level>, --log <log_level>: Verbosity levels, (default: INFO).


# Edit the reference file
python $YAG/scripts/update_assembly_edits_and_breaks.py \
    -i $FASTA \
    -o $NAME.closed.fa \
    -e $OUTDIR/$NAME.edits.txt
echo "closed gap"