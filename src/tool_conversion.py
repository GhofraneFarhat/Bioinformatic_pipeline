import argparse

def get_conversion_class(tool_name, tool_version):
    # Import the appropriate class based on the tool_name
    if tool_name == 'classify':
        from tool_conversion.classify_conversion import FastaToCsv
        return FastaToCsv
    elif tool_name == 'plASgraph':
        from tool_conversion.plasgraph_conversion import CsvToCsv
        return CsvToCsv
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

