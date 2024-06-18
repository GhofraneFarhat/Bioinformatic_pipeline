from .plasbin_flow_utils import generate_content_gc_file
from .plasbin_flow_utils import gzip_file
from .plaspipe_utils import csv_to_tsv

import os

def get_command(method_config, input_file, output_file, plasmid_scores_file = ""):

    # Extract tool_name, tool_version, script_path, tool_parameter
    tool_name = method_config['name']
    tool_version = method_config['version']
    script_name = method_config['path_to_script']
    tool_parameters = method_config['parameters']

    print(f'the plasmid score file is {plasmid_scores_file}')


    # Check if the tool is supported and return the correct command
    #plasgraph 1.0.0
    if tool_name == "plASgraph" and tool_version == "1.0.0":
        command = f'python {script_name} -i {input_file} -o {output_file}'

    #plasgraph2.0.0
    elif tool_name == "plASgraph2" and tool_version == "2.0.0":
        command = f'python {script_name} gfa {input_file} {tool_parameters} {output_file}'

    #platon1.0.0
    elif tool_name == "platon" and tool_version == "1.0.0":
        command = f'platon {input_file} {output_file}'

    #classify1.0.0
    elif tool_name == "classify" and tool_version == "1.0.0":
        command = f'python {script_name} -i {input_file} -o {output_file}'

    #bin_tool1.0.0
    elif tool_name == "bin_tool" and tool_version == "1.0.0":
        command = f'python {script_name} {input_file} {output_file}'

    #plasbin_flow1.0.0
    elif tool_name == "plasbin_flow" and tool_version == "1.0.0":
        
        #generate all the tool parameter and files
        
        #extract the outdir path and output file
        path_to_outdir, output_file = os.path.split(output_file)

        #generate the contig content gc file
        generate_content_gc_file(method_config, input_file)

        #need to add a function to verify if the file is full or not
        gc_content_file = os.path.join(method_config['plasbin_out_dir'], method_config['sample_name'] + ".gc.tsv")
        print(f'this is my gc_content_file {gc_content_file}')

        #format the input file to a gzzeped file, plasbin_flow needs it
        input_file = gzip_file(input_file)
        print(f'this is my input_file gzipped {input_file}')

        #formatting the csv plasmid score file to the plasbin-flow input (tsv)
        classification_score_result = csv_to_tsv(plasmid_scores_file)
        print(f'this is my plasmid_score_file for the plasbin_tool {classification_score_result}')


        #log_file from the yaml file
        log_file = method_config['log_file']

        command = f'python {script_name} -ag {input_file} -gc {gc_content_file} -score {classification_score_result} -out_dir {path_to_outdir} -out_file {output_file} -log_file {log_file}'


    
    else:
        raise ValueError(f"Unsupported tool: {tool_name} or version: {tool_version}")



    print(f'this is our classification file result in csv format {plasmid_scores_file}')

    return command

# Example usage
if __name__ == "__main__":

    command = get_command(config, input_file, output_file)
    print(command)