#!/usr/bin/env python
"""
This module provides functionality to dynamically load and run conversion classes
for various bioinformatics tools. It includes functions to get the appropriate
conversion class based on the tool name and version, and to run the conversion
process.

Functions:
    get_conversion_class(tool_name, tool_version): Retrieves the conversion class for a specified tool.
    run_conversion(tool_name, tool_version, input_file, output_file): Runs the conversion process for a specified tool.

The module supports conversion for the following tools:
    - plASgraph2 (version 1.0.0)
    - plasbin_flow (version 1.0.2.2)
    - PlasForest (version 1.4.0)
    - gplas2 (version 1.1.0)
    - mlplasmid (version 2.2.0)
"""
import importlib
import sys
import os
import logging
from pathlib import Path
from .plaspipe_utils import process_error, process_exception, check_file


def get_conversion_class(tool_name, tool_version):
    """
    get the conversion class for the specified tool.

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool
    Returns:
    conversion_class
    """
    #check the tool_name and version
    if not isinstance(tool_name, str) or not isinstance(tool_version, str):
        raise ValueError("Tool name and version must be strings")

    conversion_classes = {
        ('plASgraph2', '1.0.0'): ('pipeline.plaspipe.src.tools_conversion.plasgraph2_conversion', 'Plasgraph2ToCsv'),
        ('plasbin_flow', '1.0.2.2'): ('pipeline.plaspipe.src.tools_conversion.plasbin_flow_conversion', 'TsvToCsv'),
        ('PlasForest', '1.4.0'): ('pipeline.plaspipe.src.tools_conversion.PlasForest_conversion', 'PlasForestToCsv'),
        ('gplas2', '1.1.0'): ('pipeline.plaspipe.src.tools_conversion.gplas2_conversion', 'gplasToCsv'),
        ('mlplasmid', '2.2.0'): ('pipeline.plaspipe.src.tools_conversion.mlplasmid_conversion', 'mlplasmidToCsv'),
    }

    try:
        module_name, class_name = conversion_classes.get((tool_name, tool_version), (None, None))
        if module_name is None:
            raise ValueError(f"unsupported tool: {tool_name} or unsupported version: {tool_version}")

        module = importlib.import_module(module_name)
        return getattr(module, class_name)

    except ImportError as e:
        print(f"Error importing conversion class for {tool_name}: {e}")
        raise
    except AttributeError as e:
        print(f"conversion class not found for {tool_name}: {e}")
        raise
    except Exception as e:
        print(f"unexpected error in get_conversion_class: {e}")
        raise

def run_conversion(tool_name, tool_version, input_file, output_file):
    """
    Run the conversion process for the specified tool.

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool
    input_file (str): Path to the input file
    output_file (str): Path to the output file
    Returns:
    result (str): Path to the converted file
    """

    if not all(isinstance(arg, str) for arg in [tool_name, tool_version, input_file, output_file]):
        raise ValueError("All arguments must be strings")

    try:
        check_file(input_file)  # check the input file 
        
        # get the correct conversion class
        ConversionClass = get_conversion_class(tool_name, tool_version)

        conversion = ConversionClass(input_file, output_file)
        
        result = conversion.convert()
        
        check_file(result)  # check the result file 
        logging.info(f"conversion done, Output file: {result}")
        return result

    except FileNotFoundError as e:
        process_error(f"input file not found: {e}")
        raise
    except ValueError as e:
        process_error(f"conversion error: {e}")
        raise
    except Exception as e:
        process_exception(f"unexpected error in run_conversion: {e}")
        raise