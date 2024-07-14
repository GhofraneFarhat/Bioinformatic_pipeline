#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 /path/to/your/inputfile.fasta /path/to/your/outputfile.csv"
    exit 1
fi

# Assign arguments to variables
INPUT_FILE=$1
OUTPUT_FILE=$2

echo "Step 1: Checking conda installation"
which conda
if [ $? -ne 0 ]; then
    echo "Conda is not installed or not in PATH"
    exit 1
fi

echo "Step 2: Checking conda version"
conda --version

echo "Step 3: Attempting to create conda environment"
conda create -n plastest python=3.8 -y
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
conda activate plastest
if [ $? -ne 0 ]; then
    echo "Failed to activate conda environment"
    exit 1
fi

echo "Step 6: Checking Python version in the environment"
python --version

echo "Step 7: Listing installed packages in the environment"
pip list

echo "Environment setup completed successfully"

#install requirement

pip install 'scikit-learn==0.22.2.post1'
pip install pandas
pip install 'biopython==1.81'
pip install joblib
pip install requests
pip install 'numpy==1.23.0'
pip install 'Cython==0.29.14'

echo "Step 7: Checking Python requirements installed"



# Change directory to PlasForest
echo "Changing directory to PlasForest..."
cd submodules/PlasForest/
echo "Current directory: $(pwd)"

# Extract the compressed file
echo "Extracting plasforest.sav.tar.gz..."
tar -zxvf plasforest.sav.tar.gz
echo "Extraction complete."

# Make the database downloader script executable
echo "Making database downloader script executable..."
chmod 755 database_downloader.sh
echo "Permissions updated."

# Run the database downloader script
echo "Downloading PlasForest database..."
./database_downloader.sh
echo "Database download complete."

# Run PlasForest.py with the provided input and output paths
echo "Running PlasForest..."
python3 PlasForest.py -i "$INPUT_FILE" -o "$OUTPUT_FILE"
echo "PlasForest execution complete."

# Deactivate the conda environment
echo "Deactivating conda environment..."
conda deactivate
echo "Conda environment deactivated."
