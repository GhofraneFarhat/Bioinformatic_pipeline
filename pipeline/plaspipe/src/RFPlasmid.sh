#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 /path/to/your/input/folder /path/to/your/output/folder absolute/path"
    exit 1
fi

# Assign arguments to variables
INPUT_FOLDER=$1
OUTPUT_FOLDER=$2
ABSOLUTE_PATH=$3

echo "Step 1: Checking conda installation"
which conda
if [ $? -ne 0 ]; then
    echo "Conda is not installed or not in PATH"
    exit 1
fi

echo "Step 2: Checking conda version"
conda --version

echo "Step 3: Attempting to create conda environment"
conda create -n rfplasmid python=3.8 -y
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment"
    exit 1
fi

echo "Step 4: Sourcing conda.sh"
source ~/miniconda3/etc/profile.d/conda.sh
if [ $? -ne 0 ]; then
    echo "Failed to source conda.sh"
    exit 1
fi

echo "Step 5: Activating conda environment"
conda activate rfplasmid
if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment"
    exit 1
fi

echo "Step 6: Checking Python version in the environment"
python --version

echo "Step 7: Listing installed packages in the environment"
pip list

echo "Environment setup completed successfully"

# Install required packages
pip install pandas numpy scipy pysam checkm-genome biopython

# Create a folder for rfplasmid
mkdir -p ${ABSOLUTE_PATH}/submodules/plasmidRF/checkm_data

# Change directory to rfplasmid
echo "Changing directory to RFPlasmid..."
cd ${ABSOLUTE_PATH}/submodules/plasmidRF
echo "Current directory: $(pwd)"

# Change directory to checkm_data
echo "Changing directory to checkm_data..."
cd checkm_data
echo "Current directory: $(pwd)"

# Download and extract CheckM data
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar xzvf checkm_data_2015_01_16.tar.gz

# Set CheckM data root
checkm data setRoot ${ABSOLUTE_PATH}/submodules/plasmidRF/checkm_data

# Change directory to rfplasmid
cd ${ABSOLUTE_PATH}/submodules/plasmidRF
echo "Current directory: $(pwd)"


# Change directory to RFPlasmid and run setup
cd RFPlasmid/
bash getdb.sh

# Run rfplasmid.py with the provided input and output paths
echo "Running RFPlasmid..."
python3 rfplasmid.py --species Campylobacter --input "$INPUT_FOLDER" --jelly --threads 8 --out "$OUTPUT_FOLDER"
echo "RFPlasmid execution complete."

# Deactivate the conda environment
echo "Deactivating conda environment..."
conda deactivate
echo "Conda environment deactivated."