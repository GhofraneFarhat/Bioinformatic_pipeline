import argparse

def get_conversion_class(tool_name, tool_version):
    # Import the appropriate class based on the tool_name
    if tool_name == 'classify':
        from tool_conversion.classify_conversion import FastaToCsv
        return FastaToCsv
    elif tool_name == 'plASgraph':
        from tool_conversion.plasgraph_conversion import CsvToCsv
        return CsvToCsv

    elif tool_name == 'plASgraph2':
        from tool_conversion.plasgraph2_conversion import Plasgraph2ToCsv
        return Plasgraph2ToCsv

    elif tool_name == 'bin_tool':
        from tool_conversion.plasgraph_conversion import CsvToCsv
        return CsvToCsv

    elif tool_name == 'plasbin_flow':
        from tool_conversion.plasbin_flow_conversion import TsvToCsv
        return TsvToCsv

    else:
        raise ValueError(f"Unsupported tool: {tool_name}")

def run_conversion(tool_name, tool_version, input_file, output_file):
    # Get the appropriate conversion class
    ConversionClass = get_conversion_class(tool_name, tool_version)

    # Instantiate and run the conversion
    conversion = ConversionClass(input_file, output_file)
    resultat = conversion.convert()
    print(f"Conversion completed. Output file: {resultat}")
    return resultat

