import argparse
from .plaspipe_utils import process_error
from .plaspipe_utils import process_exception

def get_conversion_class(tool_name, tool_version):
    """
    Get the appropriate conversion class based on the tool name and version

    Args:
    tool_name (str): Name of the tool
    tool_version (str): Version of the tool

    Returns:
    class: The appropriate conversion class

    """
    try:
        if tool_name == 'classify' and tool_version == "1.0.0":
            from tool_conversion.classify_conversion import FastaToCsv
            return FastaToCsv

        elif tool_name == 'plASgraph' and tool_version == "1.0.0":
            from tool_conversion.plasgraph_conversion import CsvToCsv
            return CsvToCsv

        elif tool_name == 'plASgraph2' and tool_version == "2.0.0":
            from tool_conversion.plasgraph2_conversion import Plasgraph2ToCsv
            return Plasgraph2ToCsv

        elif tool_name == 'bin_tool' and tool_version == "1.0.0":
            from tool_conversion.plasgraph_conversion import CsvToCsv
            return CsvToCsv

        elif tool_name == 'plasbin_flow' and tool_version == "1.0.0":
            from tool_conversion.plasbin_flow_conversion import TsvToCsv
            return TsvToCsv

        else:
            raise ValueError(f"Unsupported tool: {tool_name} or Unsupported version: {tool_version}")

    except ImportError as e:
        process_error(f"Error importing conversion class for {tool_name}: {e}")

    except Exception as e:
        process_exception(f"Unexpected error in get_conversion_class: {e}")

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

    """
    try:
        # Get the appropriate conversion class
        ConversionClass = get_conversion_class(tool_name, tool_version)

        # Instantiate and run the conversion
        conversion = ConversionClass(input_file, output_file)
        result = conversion.convert()
        print(f"Conversion completed. Output file: {result}")
        return result

    except ValueError as e:
        process_error(f"Conversion error: {e}")
    except Exception as e:
        process_exception(f"Unexpected error in run_conversion: {e}")

