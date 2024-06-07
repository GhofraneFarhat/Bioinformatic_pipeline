def get_command(method_config, input_file, output_file):
    # Define command templates for each tool
    command_tool = {
        'plASgraph': 'python {script} -i {input_file} -o {output_file}',
        'plagraph2': 'python {script} --input {input_file} --output {output_file}',
        'platon': 'platon {input_file} {output_file}',
        'classify': 'python {script} -i {input_file} -o {output_file}'
    }

    # Extract tool_name, tool_version, script_path
    tool_name = method_config['name']
    tool_version = method_config['version']
    script_name = method_config['command']
    tool_parameters = method_config['parameters']

    # Check if the tool is supported
    if tool_name not in command_tool:
        raise ValueError(f"Unsupported tool: {tool_name}")

    # Get the command template for the specified tool
    command_template = command_tool[tool_name]

    # Format the command with the provided parameters
    command = command_template.format(script=script_name, input_file=input_file, output_file=output_file)

    # Add additional parameters if provided


    return command

# Example usage
if __name__ == "__main__":

    command = get_command(config, input_file, output_file)
    print(command)