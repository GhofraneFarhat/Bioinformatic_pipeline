import argparse
import logging
from importlib import import_module
from .plaspipe_utils import process_error, process_exception, check_file

def get_conversion_class(tool_name, tool_version):
    """
    Get the appropriate conversion class based on the tool name and version.

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool

    Returns:
    class: The appropriate conversion class
    """
    
    if not isinstance(tool_name, str) or not isinstance(tool_version, str):
        raise ValueError("Tool name and version must be strings")

    conversion_classes = {
        ('classify', '1.0.0'): ('pipeline.plaspipe.tool_conversion.classify_conversion', 'FastaToCsv'),
        ('plASgraph', '1.0.0'): ('pipeline.plaspipe.tool_conversion.plasgraph_conversion', 'CsvToCsv'),
        ('plASgraph2', '2.0.0'): ('pipeline.plaspipe.tool_conversion.plasgraph2_conversion', 'Plasgraph2ToCsv'),
        ('bin_tool', '1.0.0'): ('pipeline.plaspipe.tool_conversion.plasgraph_conversion', 'CsvToCsv'),
        ('PlasBin', '1.0.0'): ('pipeline.plaspipe.tool_conversion.PlasBin_conversion', 'PlasBinToCsv'),
        ('plasbin_flow', '1.0.0'): ('pipeline.plaspipe.tool_conversion.plasbin_flow_conversion', 'TsvToCsv'),
        ('PlasForest', '1.4.0'): ('pipeline.plaspipe.tool_conversion.PlasForest_conversion', 'PlasForestToCsv'),
        ('gplas2', '1.1.0'): ('pipeline.plaspipe.tool_conversion.gplas2_conversion', 'gplasToCsv'),
    }

    try:
        module_name, class_name = conversion_classes.get((tool_name, tool_version), (None, None))
        if module_name is None:
            raise ValueError(f"unsupported tool: {tool_name} or unsupported version: {tool_version}")

        module = import_module(module_name)
        return getattr(module, class_name)

    except ImportError as e:
        process_error(f"error importing conversion class for {tool_name}: {e}")
        raise
    except AttributeError as e:
        process_error(f"conversion class not found for {tool_name}: {e}")
        raise
    except Exception as e:
        process_exception(f"unexpected error in get_conversion_class: {e}")
        raise

def run_conversion(tool_name, tool_version, input_file, output_file):
    """
    run the conversion process for the specified tool.

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