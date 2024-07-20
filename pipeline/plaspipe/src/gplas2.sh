#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 /path/to/the/inputfile.gfa path/to/classification/file/tab outputfile/name min_length"
    exit 1
fi

# Assign arguments to variables
INPUT_FILE="$1"
CLASSIFICATION_RESULT="$2"
OUTPUT_NAME="$3"
MIN_LENGTH="$4"

# Change directory to gplas (gplas is a submodule in submodules folder)
echo "Changing directory to gplas2..."
cd submodules/gplas2/ || { echo "Failed to change directory to gplas2"; exit 1; }
echo "Current directory: $(pwd)"

echo "Step 1: Checking conda installation"
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or not in PATH"
    exit 1
fi

echo "Step 2: Checking conda version"
conda --version

echo "Step 3: Attempting to create conda environment"
if ! conda env create --name gplas2 --file envs/gplas.yaml -y; then
    echo "Failed to create conda environment"
    exit 1
fi

echo "Step 4: Sourcing conda.sh"
CONDA_SH_PATH="$HOME/miniconda3/etc/profile.d/conda.sh"
if [ -f "$CONDA_SH_PATH" ]; then
    source "$CONDA_SH_PATH"
else
    echo "conda.sh not found at $CONDA_SH_PATH"
    exit 1
fi

echo "Step 5: Activating conda environment"
if ! conda activate gplas2; then
    echo "Failed to activate conda environment"
    exit 1
fi

echo "Step 6: Running setup for gplas2"
if ! pip install -e .; then
    echo "Failed to install gplas2"
    exit 1
fi

echo "Successfully installed"

# Run gplas2 with the provided input and output paths
echo "Running gplas2..."
if ! gplas -i "$INPUT_FILE" -c extract -n "$OUTPUT_NAME" -l "$MIN_LENGTH"; then
    echo "gplas2 extraction failed"
    exit 1
fi
echo "Running gplas2 extraction complete."

if ! gplas -c predict -i "$INPUT_FILE" -P "$CLASSIFICATION_RESULT" -n "$OUTPUT_NAME"; then
    echo "gplas2 prediction failed"
    exit 1
fi
echo "gplas2 execution complete."

# Deactivate the conda environment
echo "Deactivating conda environment..."
conda deactivate
echo "Conda environment deactivated."