from .plasbin_flow_utils import generate_content_gc_file
from .plasbin_utils import generate_seed_contigs
from .plasbin_utils import run_blast6
from .plasbin_flow_utils import gzip_file

from .plaspipe_utils import csv_to_tsv
from .plaspipe_utils import check_file
from .plaspipe_utils import process_exception
from .plaspipe_utils import process_error
from .plaspipe_utils import process_arguments
from .plaspipe_utils import absolute_path
from .plaspipe_utils import conversion_gfa_fasta
from .plaspipe_utils import import_file
from .plaspipe_utils import csv_to_tab
from .plaspipe_utils import update_contig_names
from .plaspipe_utils import gzip_file


import os
import logging
import subprocess
import sys

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
        python_executable = sys.executable

        # Generate command based on tool and version

        if tool_name == "plASgraph2" and tool_version == "2.0.0":
            #the script_path
            plasgraph2_script = os.path.join(absolute_path(), 'submodules/plASgraph2/src/plASgraph2_classify.py')
            # Use sys.executable to get the correct Python interpreter
             
            
            command = [
                python_executable,
                plasgraph2_script,
                'gfa',
                input_file,
                tool_parameters,
                output_file
            ]


        elif tool_name == "platon" and tool_version == "1.0.0":
            command = f'platon {input_file} {output_file}'


        elif tool_name == "plasbin_flow" and tool_version == "1.0.0":
            # Handle plasbin_flow specific requirements
            path_to_outdir, output_file = os.path.split(output_file)

            sample_name = method_config['sample_name']
            log_file = method_config.get('log_file', 'plasbin_flow.log')  # Provide a default log file name
            assembler = method_config['assembler']
            seed_score = method_config['seed_score']
            seed_len = method_config['seed_len']
            gc_intervals = method_config['gc_interval_file']
            out_gc_file = method_config['plasbin_out_dir']

            #the default values
            default_gc_intervals = os.path.join(absolute_path(), 'submodules/PlasBin-flow/example/default/gc_intervals.txt')
            default_out_dir = os.path.join(absolute_path(), 'out/plasbin_flow')

            #the script_path
            plasbin_flow_script = os.path.join(absolute_path(), 'submodules/PlasBin-flow/code/plasbin_flow.py')
            
            #processing the attributes
            argument = [sample_name,assembler, seed_score, seed_len, gc_intervals, out_gc_file]
            default_arg = ['test1','unicycler', 0.58, 2650, default_gc_intervals,default_out_dir]
    
            plasbin_flow_argument = process_arguments(argument, default_arg)
            sample_name, assembler, seed_score, seed_len, gc_intervals, out_gc_file = plasbin_flow_argument



            # Generate GC content file
            try:
                generate_content_gc_file(method_config, input_file)
                gc_content_file = os.path.join(out_gc_file, sample_name + ".gc.tsv")
                
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

 
              
            #get the command
            command = [
                python_executable,
                plasbin_flow_script,
                '-ag', input_file,
                '-gc', gc_content_file,
                '-score', classification_score_result,
                '-out_dir', path_to_outdir,
                '-out_file', output_file,
                '-log_file', log_file,
                '-seed_len', str(seed_len),
                '-seed_score', str(seed_score),
                '-assembler', assembler,
                '-gc_intervals', gc_intervals
            ]


        elif tool_name == "PlasForest" and tool_version == "1.4.0":

            # Define the path to the bash script
            bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/PlasForest.sh')

            # Make sure the bash script is executable
            os.chmod(bash_script, 0o755)

            # create the command to run the bash script with the input and output file paths as arguments
            command = [bash_script, input_file, output_file]


        elif tool_name == "RFPlasmid" and tool_version == "1.0.0":

            project_path = absolute_path()
            path_to_outdir, output_file = os.path.split(output_file)
            # Define the path to the bash script
            bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/RFPlasmid.sh')


            # Create the folder for the input files for RFPlasmid
            fasta_folder = os.path.join(project_path, "fasta_folder")
            os.makedirs(fasta_folder, exist_ok=True)

            # Import the input file to the fasta folder
            import_file(input_file, fasta_folder)

            # Make sure the bash script is executable
            os.chmod(bash_script, 0o755)

            # create the command to run the bash script with the input and output file paths as arguments
            command = [bash_script, fasta_folder, path_to_outdir, project_path]


        elif tool_name == "gplas2" and tool_version == "1.1.0":

            # Define the path to the bash script
            bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/gplas2.sh')
            min_length = method_config['min_length']
            output_name = method_config['output_name']   
  

            #default parameter
            argument = [min_length, output_name]
            default_arg = [100,'gplas']
    
            gplas_argument = process_arguments(argument, default_arg)
            min_length, output_name = gplas_argument

            # Make sure the bash script is executable
            os.chmod(bash_script, 0o755)

            #create a tab file for gplas2
            try:
                classification_result = csv_to_tab(plasmid_scores_file)
                check_file(classification_result)  # Verify the converted file exists
            except Exception as e:
                process_error(f"Error converting plasmid scores file to TAB: {e}")

            #change contigs name for gplas2
            try:
                classification_result_file = update_contig_names(classification_result, input_file)
                check_file(classification_result_file)  # Verify the converted file exists
            except Exception as e:
                process_error(f"Error updating contigs name: {e}")


            # create the command to run the bash script with the input and output file paths as arguments
            command = [
                "bash", 
                bash_script, 
                input_file, 
                classification_result_file, 
                output_name, 
                str(min_length)
                ]


        elif tool_name == "mlplasmid" and tool_version == "2.2.0":

            path_to_outdir, output_file_name = os.path.split(output_file)

            #define default values
            threshold = method_config['threshold']
            species = method_config['species']

            #check species
            mlplasmid_species = ['Enterococcus faecium', 'Klebsiella pneumoniae', 'Escherichia coli']
            
            if species not in mlplasmid_species:
                process_error(f"Unsupported species '{species}'. Supported species are: {', '.join(mlplasmid_species)}")

            #set the default threshold
            if threshold == None:
                threshold = 0.7

            #gzip the input file 
            input_file = gzip_file(input_file, path_to_outdir)

            # Define the path to the bash script
            bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/mlplasmid.sh')

            # Make sure the bash script is executable
            os.chmod(bash_script, 0o755)

            # create the command to run the bash script with the input and output file paths as arguments
            command = [bash_script, input_file, output_file, str(threshold), species]      


        else:
            raise ValueError(f"Unsupported tool: {tool_name} or version: {tool_version}")

        return command

    except KeyError as e:
        process_error(f"Missing key in method_config: {e}")
    except Exception as e:
        process_exception(f"Unexpected error in get_command: {e}")

