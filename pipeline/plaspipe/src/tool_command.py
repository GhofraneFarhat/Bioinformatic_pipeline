from .plasbin_flow_utils import generate_content_gc_file
from .Plasbin_utils import generate_seed_contigs
from .plasbin_flow_utils import gzip_file

from .plaspipe_utils import csv_to_tsv
from .plaspipe_utils import check_file
from .plaspipe_utils import process_exception
from .plaspipe_utils import process_error
from .plaspipe_utils import process_arguments
from .plaspipe_utils import absolute_path


import os
import logging

def get_command(method_config, input_file, output_file, plasmid_scores_file=""):
    """
    Generate the command for running various bioinformatics tools based on the method configuration.

    Args:
    method_config (dict): Configuration dictionary for the method.
    input_file (str): Path to the input file.
    output_file (str): Path to the output file.
    plasmid_scores_file (str): Path to the plasmid scores file (optional).

    Returns:
    str: The generated command to run the tool.

    Raises:
    ValueError: If the tool or version is not supported.
    """
    try:
        # Extract tool information from method_config
        tool_name = method_config['name']
        tool_version = method_config['version']
        tool_parameters = method_config['parameters']

        logging.info(f"Generating command for {tool_name} version {tool_version}")
        logging.debug(f"Input file: {input_file}")
        logging.debug(f"Output file: {output_file}")

        # Generate command based on tool and version
        if tool_name == "plASgraph" and tool_version == "1.0.0":

            #the script_path
            plasgraph_script = os.path.join(absolute_path(), 'submodules\plASgraph\plASgraph.py')
            
            command = f'python {script_name} -i {input_file} -o {output_file}'


        elif tool_name == "plASgraph2" and tool_version == "2.0.0":

            #the script_path
            plasgraph2_script = os.path.join(absolute_path(), 'submodules\plASgraph2\src\plASgraph2_classify.py')
            
            command = f'python {plasgraph2_script} gfa {input_file} {tool_parameters} {output_file}'


        elif tool_name == "platon" and tool_version == "1.0.0":
            command = f'platon {input_file} {output_file}'


        elif tool_name == "classify" and tool_version == "1.0.0":

            classify_script_name = os.path.join(absolute_path(), 'submodules\classify\classify.py')
            
            command = f'python {classify_script_name} -i {input_file} -o {output_file}'


        elif tool_name == "bin_tool" and tool_version == "1.0.0":
            classify_script_name = os.path.join(absolute_path(), 'submodules\bin_tool\bin.py')
            command = f'python {script_name} {input_file} {output_file}'


        elif tool_name == "plasbin_flow" and tool_version == "1.0.0":
            # Handle plasbin_flow specific requirements
            path_to_outdir, output_file = os.path.split(output_file)

            # Generate GC content file
            try:
                generate_content_gc_file(method_config, input_file)
                gc_content_file = os.path.join(method_config['plasbin_out_dir'], method_config['sample_name'] + ".gc.tsv")
                check_file(gc_content_file)  # Verify the file was created successfully
                logging.info(f"GC content file is saved to {gc_content_file}")
            except Exception as e:
                process_error(f"Error generating or accessing GC content file: {e}")


            # Gzip input file
            try:
                input_file = gzip_file(input_file)
                check_file(input_file)  # Verify the gzipped file exists
                logging.info(f"Gzziped file {input_file}")
            except Exception as e:
                process_error(f"Error gzipping input file: {e}")


            # Convert CSV plasmid score file to TSV
            try:
                classification_score_result = csv_to_tsv(plasmid_scores_file)
                check_file(classification_score_result)  # Verify the converted file exists
            except Exception as e:
                process_error(f"Error converting plasmid scores file to TSV: {e}")


            log_file = method_config.get('log_file', 'plasbin_flow.log')  # Provide a default log file name
            assembler = method_config['assembler']
            seed_score = method_config['seed_score']
            seed_len = method_config['seed_len']
            gc_intervals = method_config['gc_interval_file']

            #the default values
            default_gc_intervals = os.path.join(absolute_path(), 'submodules\PlasBin-flow\example\default\gc_intervals.txt')
            
            #the script_path
            plasbin_flow_script = os.path.join(absolute_path(), 'submodules\PlasBin-flow\code\plasbin_flow.py')
            
            #processing the attributes
            argument = [assembler, seed_score, seed_len, gc_intervals]
            default_arg = ['unicycler', 0.58, 2650, default_gc_intervals ]
    
            plasbin_flow_argument = process_arguments(argument, default_arg)
            assembler, seed_score, seed_len, gc_intervals = plasbin_flow_argument 
              
            #get the command
            command = f'python {plasbin_flow_script} -ag {input_file} -gc {gc_content_file} -score {classification_score_result} -out_dir {path_to_outdir} -out_file {output_file} -log_file {log_file} -seed_len {seed_len} -seed_score {seed_score} -assembler {assembler} -gc_intervals {gc_intervals}'
        
        elif tool_name == "PlasBin" and tool_version == "1.0.0":

            #get the output folder
            path_to_outdir, output_file = os.path.split(output_file)

            mapping_file = method_config['mapping_file']
            alpha1 = method_config['alpha1']
            alpha2 = method_config['alpha2']
            rmiter = method_config['rmiter']

            #assign default value
            default_maping_file = os.path.join(absolute_path(), 'submodules/PlasBin/example/test/input/filtered_genes_to_contigs.csv')
            
            argument = [mapping_file, alpha1, alpha2, rmiter]
            default_arg = [default_maping_file, 1, 1, 50]
    
            plasbin_argument = process_arguments(argument, default_arg)
            mapping_file, alpha1, alpha2, rmiter = plasbin_argument

            try:
                generate_seed_contigs(input_file, mapping_file, method_config)
                seed_contigs_file = os.path.join(method_config['plasbin_out_dir'], "seed_contigs.csv")
                check_file(seed_contigs_file)  # Verify the file was created successfully
                logging.info(f"seed contig file is saved to {seed_contigs_file}")
            except Exception as e:
                process_error(f"Error generating or accessing seed contigs file: {e}")

            #the script_path
            plasbin_script = os.path.join(absolute_path(), 'submodules\PlasBin\code\plasmids_iterative.py')
            
            command = f'python {plasbin_script} --ag {input_file} --map {mapping_file} --seeds {seed_contigs_file} --out {path_to_outdir}'


        else:
            raise ValueError(f"Unsupported tool: {tool_name} or version: {tool_version}")

        return command

    except KeyError as e:
        process_error(f"Missing key in method_config: {e}")
    except Exception as e:
        process_exception(f"Unexpected error in get_command: {e}")

