# Tool Conversion Module

## Overview

The `tools_conversion` module is part of the Bioinformatic Pipeline project. It contains classes and functions designed to convert output files from various bioinformatics tools into a standardized CSV format (pipeline format). This module simplifies the integration of different tool outputs into a unified analysis pipeline.

## Structure

The `tools_conversion` directory contains the following files:

- `__init__.py`: Initializes the `tools_conversion` package.
- `plasgraph2_conversion.py`: Contains the `Plasgraph2ToCsv` class for converting Plasgraph2 output file.
- `plasbin_flow_conversion.py`: Contains the `PlasbinFlowToCsv` class for converting Plasbin Flow output file.
- `PlasForest_conversion.py`: Contains the `PlasForestToCsv` class for converting PlasForest output file.
- `gplas2_conversion.py`: Contains the `gplasToCsv` class for converting gplas2 output file.
- `mlplasmid_conversion.py`: Contains the `mlplasmidToCsv` class for converting mlplasmid output file.

## Usage

To use the classes in this module, you can import them in your script. Here’s an example of how to use the `Plasgraph2ToCsv` class:

```python
from pipeline.plaspipe.src.tools_conversion.plasgraph2_conversion import Plasgraph2ToCsv

# Define input and output file
input_file = 'path/to/input/plasgraph2_output'
output_file = 'path/to/output/plasgraph2_converted_output.csv'

# Create an instance of the conversion class
converter = Plasgraph2ToCsv(input_file, output_file)

# Conversion
converter.convert()