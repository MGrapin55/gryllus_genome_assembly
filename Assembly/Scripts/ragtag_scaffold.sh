#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --job-name=ragtag_scaffold
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20gb
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
#SBATCH --qos=short
set -euo pipefail
# ENV
module purge
module load anaconda/25.3

set +eou pipefail
source ~/.bashrc
conda activate "$NRDSTOR/ragtag"
set -eou pipefail

# ====================================================================================================== #
# Purpose: Homology-based assembly scaffolding: Order and orient sequences in 'query.fa' by comparing
# them to sequences in 'reference.fa'
# ====================================================================================================== #
# Directory Location 
OUTDIR=$1

# Path to reference fasta
REF="/work/moorelab/mgrapin2/GRYLLUS_ASSEMBLIES/ref/assimilis_chr_reference.fa"

# Path to query fasta
QUERY=$2

# Iteration of Scaffolding 
V=1

##########################################################################################################
##											Commands    												##
##########################################################################################################
mkdir -p $OUTDIR

# scaffold a query assembly
ragtag.py scaffold $REF $QUERY -o $OUTDIR -t $SLURM_CPUS_PER_TASK -r -u
# (optional parameters to play with if you want more strict)  -i 0.95 -a 0.95 -s 0.95

# scaffolding options:
#  -e <exclude.txt>     list of reference sequences to ignore [null]
#  -j <skip.txt>        list of query sequences to leave unplaced [null]
#  -J <hard-skip.txt>   list of query headers to leave unplaced and exclude from 'chr0' ('-C') [null]
#  -f INT               minimum unique alignment length [1000]
#  --remove-small       remove unique alignments shorter than '-f'
#  -q INT               minimum mapq (NA for Nucmer alignments) [10]
#  -d INT               maximum alignment merge distance [100000]
#  -i FLOAT             minimum grouping confidence score [0.2]
#  -a FLOAT             minimum location confidence score [0.0]
#  -s FLOAT             minimum orientation confidence score [0.0]
#  -C                   concatenate unplaced contigs and make 'chr0'
#  -r                   infer gap sizes. if not, all gaps are 100 bp
#  -g INT               minimum inferred gap size [100]
#  -m INT               maximum inferred gap size [100000]

# input/output options:
#  -o PATH              output directory [./ragtag_output]
#  -w                   overwrite intermediate files
#  -u                   add suffix to unplaced sequence headers

# mapping options:
#  -t INT               number of minimap2/unimap threads [1]
#  --aligner PATH       aligner executable ('nucmer', 'unimap' or 'minimap2') [minimap2]
#  --mm2-params STR     space delimited minimap2 parameters (overrides '-t') ['-x asm5']
#  --unimap-params STR  space delimited unimap parameters (overrides '-t') ['-x asm5']
#  --nucmer-params STR  space delimted nucmer parameters ['--maxmatch -l 100 -c 500']

echo "######################################################################"
echo "=== SLURM Job Efficiency Report ==="
echo "Check with: seff $SLURM_JOB_ID"
echo "######################################################################"

echo "Job finished at: $(date)"
