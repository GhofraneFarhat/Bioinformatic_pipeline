import argparse
import logging
from .plaspipe_utils import process_error, process_exception, check_file

def get_conversion_class(tool_name, tool_version):

    """
    Get the appropriate conversion class based on the tool name and version

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool

    Returns:
    class: The appropriate conversion class

    Raises:
    ValueError: If the tool or version is not supported
    ImportError: If there's an error importing the conversion class
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
    }

    try:
        module_name, class_name = conversion_classes.get((tool_name, tool_version), (None, None))
        if module_name is None:
            raise ValueError(f"Unsupported tool: {tool_name} or Unsupported version: {tool_version}")

        module = __import__(module_name, fromlist=[class_name])
        return getattr(module, class_name)

    except ImportError as e:
        process_error(f"Error importing conversion class for {tool_name}: {e}")
        raise
    except AttributeError as e:
        process_error(f"Conversion class not found for {tool_name}: {e}")
        raise
    except Exception as e:
        process_exception(f"Unexpected error in get_conversion_class: {e}")
        raise

def run_conversion(tool_name, tool_version, input_file, output_file):
    """
    Run the conversion process for the specified tool

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool
    input_file (str): Path to the input file
    output_file (str): Path to the output file

    Returns:
    str: Path to the converted file

    Raises:
    ValueError: If the input parameters are invalid
    FileNotFoundError: If the input file is not found
    """
    if not all(isinstance(arg, str) for arg in [tool_name, tool_version, input_file, output_file]):
        raise ValueError("All arguments must be strings")

    try:
        check_file(input_file)  # Ensure input file exists
        
        # Get the appropriate conversion class
        ConversionClass = get_conversion_class(tool_name, tool_version)

        # Instantiate and run the conversion
        conversion = ConversionClass(input_file, output_file)
        
        result = conversion.convert()
        
        
        check_file(result)  # Ensure output file was created
        logging.info(f"Conversion completed. Output file: {result}")
        return result

    except FileNotFoundError as e:
        process_error(f"Input file not found: {e}")
        raise
    except ValueError as e:
        process_error(f"Conversion error: {e}")
        raise
    except Exception as e:
        process_exception(f"Unexpected error in run_conversion: {e}")
        raise

