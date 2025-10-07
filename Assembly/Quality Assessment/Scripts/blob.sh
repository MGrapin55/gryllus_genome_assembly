#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --job-name=BlobTools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb		
#SBATCH --mail-user=
#SBATCH --mail-type=END,FAIL
#SBATCH --error=./std.err/%x.%J.err
#SBATCH --output=./std.out/%x.%J.out
#SBATCH --partition=batch
set -eou pipefail

# ENVIRONMENT SETUP
module purge
module load biodata/1.0
module load blast/2.15
module load anaconda/25.3
module load samtools/1.20

set +eou pipefail
source ~/.bashrc
conda activate "$NRDSTOR/btk"
set -eou pipefail

# ====================================================================================================== #
# Purpose: Runs BlobToolKit Steps and Prints Viewer/API Instructions
# ====================================================================================================== #

# Working Directory 
WKDIR=""


# Name of Assembly
ASSEMBLY_NAME=""

# BLAST Results
BLAST=

# BAM file 
BAM=""

BUSCO=""

# Tax Id
ID=

# NCBI Taxonomy
TAXDUMP=""

################################################################################
#                       Ensure dataset folder exists                            #
################################################################################
cd "$WKDIR"
DATASET_DIR=$WKDIR/${ASSEMBLY_NAME}_Blob_DB
if [[ ! -d "$DATASET_DIR" ]]; then
    echo "Creating dataset directory for BlobToolKit..."
    mkdir -p "$DATASET_DIR"
    
fi

################################################################################
#                       Run BlobToolKit Pipeline Steps                           #
################################################################################

# 1) Make a blob database
if [[ ! -f "$DATASET_DIR/meta.json" ]]; then
echo "Creating BlobToolKit database..."
blobtools create \
    --fasta "${ASSEMBLY_NAME}.fa" \
    --meta config.yaml \
    --taxid "$ID" \
    --taxdump "$TAXDUMP" \
    "$DATASET_DIR"
fi

# 2) Add BLAST hits
if [[ -f "$BLAST" ]]; then
    echo "Adding BLAST hits to BlobToolKit database..."
    blobtools add \
        --hits "$BLAST" \
        --taxrule bestsumorder \
        --taxdump "$TAXDUMP" \
        "$DATASET_DIR"
fi

# 3) Add coverage
if [[ -f "$BAM" ]]; then
echo "Adding coverage to BlobToolKit database..."
blobtools add \
    --cov "$BAM" \
    --threads "$SLURM_CPUS_PER_TASK" \
    "$DATASET_DIR"
fi

if [[ -f "$BUSCO" ]]; then
echo "Adding busco to BlobToolKit database..."
blobtools add \
    --busco $BUSCO
fi

echo "Blob DB: $DATASET_DIR"
################################################################################
#                       Viewer/API Instructions                                 #
################################################################################

echo ""
echo "=================================================================="
echo "BlobToolKit pipeline steps completed!"
echo ""
echo "To start the API and Viewer with standard ports, run:"
echo ""
echo "# Start the API in background:"
echo "BTK_API_PORT=8880 BTK_PORT=8881 BTK_FILE_PATH=$WKDIR ./blobtoolkit-api &"
echo ""
echo "# Start the Viewer in background:"
echo "BTK_API_PORT=8880 BTK_PORT=8881 ./blobtoolkit-viewer &"
echo ""
echo "SSH port forwarding from your local machine:"
echo "ssh -L 8881:127.0.0.1:8881 -L 8880:127.0.0.1:8880 ${USER}@swan"
echo ""
echo "Then open your browser at:"
echo "http://localhost:8881/view/all"
echo ""
echo "Alternative: if ports 8880/8881 are in use locally, use these ports instead:"
echo ""
echo "ssh -L 9001:127.0.0.1:8881 -L 9000:127.0.0.1:8880 ${USER}@swan"
echo "Then open browser at:"
echo "http://localhost:9001/view/all"
echo "=================================================================="