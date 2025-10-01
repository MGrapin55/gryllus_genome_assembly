#!/bin/bash 

# Author: Michael Grapin -- Moore Lab Research Technician 
# Email: mgrapin2@nebraska.edu

# DESCRIPTION 
# Hifi Data Genome Assembly and Annotation Pipeline 


# Directory Structure
# 00_Raw_Data
# 01_Initial_QC
# 02_Filtered_Adaptors
# 03_Count_Kmers
#   --{k}mer_Files
#                   --{k}mer.stats
#                   --{k}mer.histo      # Hisogram for GenomeScope2.0
#                   --{k}mer_DB         # Meryl Kmer Database
#                               -- GenomeScope
#                                        -- .plots           # GenomeScope Plots
# 04_Assembly
#   --common        # Shared Hifiasm Files
#   --run{N}        # Iteration of assembly specific files
# 05_QC_Assessment
#   --quast
#   --busco
#   --merqury
#   --blobtoolkit
# 06_Scaffold
#   --run{N}
#       --QC_Assessment
#           --quast

# To be determined...

##########################################################################################################
##											Analysis Setup												##
##########################################################################################################
# Usage: bash SCRIPTS/00_Analysis_Setup.sh

# Central Location of Analysis
WKDIR=

# Raw Hifi Data Files (Absolute Path) 
BAM_1=""
# Make up to BAM_N...

# ====================================================================================================== #
echo "Starting Setup..."

cd "$WKDIR"
echo "Central Working Directory: $(pwd)"
echo "# ======================================== # "
echo "Making 00_Raw_Data..."
mkdir -p "$WKDIR/00_Raw_Data"

echo "Copying BAM File(s) to 00_Raw-Data..."
[ -e "$BAM_1" ] || echo "$(basename "$BAM_1") does not exist."
rsync -ah "$BAM_1" "$WKDIR/00_Raw_Data/"
echo "Done Copying BAM File(s)..."

echo " Make std.err and std.out Directories... "
mkdir -p std.err
mkdir -p std.out
echo "Generating Pipeline Useful Information..."
echo "# Pipeline Useful Information" > setup.log
echo "# Central Working Directory: $WKDIR" >> setup.log
echo "# std.err path: #SBATCH --error=$WKDIR/std.err/%x.%J.err" >> setup.log
echo "# std.out path: #SBATCH --output=$WKDIR/std.out/%x.%J.out" >> setup.log
echo "Check: setup.log"
echo "Done Setup!"
echo "Proceed to SCRIPTS/01_FastPLong.sh"

    
