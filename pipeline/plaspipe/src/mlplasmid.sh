#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_file> <output_file> <threshold> <species>"
    exit 1
fi


INPUT_FILE=$1
OUTPUT_FILE=$2
THRESHOLD=$3
SPECIES=$4

# change directory to mlplasmid
cd submodules/mlplasmids
echo "Current directory: $(pwd)"

# Run the R script with provided arguments
echo "mlplasmid running ..."
Rscript scripts/run_mlplasmids.R "$INPUT_FILE" "$OUTPUT_FILE" "$THRESHOLD" "$SPECIES"
echo "running completed"