def get_command(method_config, input_file, output_file):
    # Define command templates for each tool
    command_tool = {
        'plASgraph': 'python {script} -i {input_file} -o {output_file}',
        'plASgraph2': 'python {script} gfa {input_file} {parameter} {output_file}',
        'platon': 'platon {input_file} {output_file}',
        'classify': 'python {script} -i {input_file} -o {output_file}',
        'bin_tool': 'python {script} {input_file} {output_file}'
    }

    # Extract tool_name, tool_version, script_path
    tool_name = method_config['name']
    tool_version = method_config['version']
    script_name = method_config['path_to_script']
    tool_parameters = method_config['parameters']

    # Check if the tool is supported
    if tool_name not in command_tool:
        raise ValueError(f"Unsupported tool: {tool_name}")

    # Get the command template for the specified tool
    command_template = command_tool[tool_name]

    # Format the command with the provided parameters
    command = command_template.format(script=script_name, input_file=input_file, parameter = tool_parameters, output_file=output_file)

    # Add additional parameters if provided


    return command

# Example usage
if __name__ == "__main__":

    command = get_command(config, input_file, output_file)
    print(command)