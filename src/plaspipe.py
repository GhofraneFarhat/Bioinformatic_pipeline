import yaml
import os
import json
import time
import argparse
import gzip
import sys
import shutil
import logging
from Bio import SeqIO

from .plaspipe_utils import process_exception
from .plaspipe_utils import check_file
from .plaspipe_utils import gunzip_FASTA
from .plaspipe_utils import gunzip_GFA
from .plaspipe_utils import create_directory


from .plaspipe_data import PipelineData



def parse_arguments():

    """
    Read arguments from the command line.
    Returns: args (list): List of arguments.
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(description='Pipeline for bioinformatic tools (classification and binning tools for plasmid analysis')
    parser.add_argument('user_file', type=str, help='Path to the YAML file provided by the user')
    parser.add_argument('pipeline_output', type=str, help='File to save contigs dict and bins dict (pipeline output)')
    args = parser.parse_args()

    # Verify the number of arguments
    if len(sys.argv) != 3:
        process_exception("Invalid number of arguments. Usage: python script.py user_file.yaml output_file.txt")

    # Verify the user_file is a YAML file
    if not args.user_file.endswith('.yaml'):
        process_exception(f"Invalid user file '{args.user_file}'. Expected a YAML file.")



    return args


# the yaml file, from the user need some changes
def load_yaml(file):

    """
    Load the YAML file.
    Args:
        file (str): Path to the YAML file.
    Returns:
        configs (dict): Configurations from the YAML file.
    """

    try:
        with open(file, 'r') as yaml_file:
            configs = yaml.safe_load(yaml_file)
        if not configs:
            process_exception(f'YAML file {file} is empty.')
        return configs
    except Exception as e:
        process_exception(f'Error loading YAML file {file}: {e}')


def get_input_file(method_config):

    outdir_pipeline = method_config['outdir_pipeline']
    os.makedirs(outdir_pipeline, exist_ok=True)

    gun_gfa_path = os.path.join(method_config['outdir_pipeline'], "gunzipped_gfa.gfa") 
    gun_fasta_path = os.path.join(method_config['outdir_pipeline'], "gunzipped_fasta.fasta")  

    gfa_path = method_config['input']['path_to_input_gfa']
    fasta_path = method_config['input']['path_to_input_fasta']
    print(f"this is my gfa file {gfa_path}")


    #unzip the gfa file
    try:
        if gfa_path and gfa_path.endswith('.gz'):
            gfa_path = gunzip_GFA(gfa_path, gun_gfa_path)
    except Exception as e:
        logging.error(f"Error unzipping GFA file {gfa_path}: {e}")


    #unzip the fasta file
    try:
        if fasta_path and fasta_path.endswith('.gz'):
            fasta_path = gunzip_FASTA(fasta_path, gun_fasta_path)
    except Exception as e:
        logging.error(f"Error unzipping FASTA file {fasta_path}: {e}")


    return gfa_path, fasta_path


def run_plaspipe_data(method_configs, prefix, class_tool_name, class_tool_version, bin_tool_name, bin_tool_version, input_gfa, input_fasta, classification_dir, binning_dir):
            
    #cretae the plaspipe_data instance
    plaspipe_data = PipelineData(method_configs, prefix, class_tool_name, class_tool_version, bin_tool_name, bin_tool_version, input_gfa, input_fasta, classification_dir, binning_dir)

    #load contigs name
    plaspipe_data.load_contigs()

    #update the plaspipe_data
    plaspipe_data.update_plaspipe_data_from_classwrapper()

    #update the plaspipe_data
    plaspipe_data.update_plaspipe_data_from_binwrapper()

    return plaspipe_data


#pipeline configs 
def get_pipeline_configs(method_configs):
    """
    Extract pipeline configurations from the method_configs dictionary.

    Args:
        method_configs (dict): Configuration dictionary for the pipeline.

    Returns:
        list: A list containing the extracted configurations.
    """
    pipeline_configs = []

    # Extract output directories
    classification_dir = method_configs['classification']['output']['outdir_classification']
    binning_dir = method_configs['binning']['output']['outdir_binning']
    prefix = method_configs['prefix']
    pipeline_configs.extend([classification_dir, binning_dir, prefix])

    # Extract tool information
    class_tool_name = method_configs['classification']['name']
    class_tool_version = method_configs['classification']['version']
    class_input_format = method_configs['classification']['input_format']
    bin_tool_name = method_configs['binning']['name']
    bin_tool_version = method_configs['binning']['version']
    bin_input_format = method_configs['binning']['input_format']
    pipeline_configs.extend([class_tool_name, class_tool_version, class_input_format,
                             bin_tool_name, bin_tool_version, bin_input_format])

    # Extract output directory for the pipeline
    out_dir = method_configs['outdir_pipeline']
    pipeline_configs.append(out_dir)

    return pipeline_configs


#the function to generate the pipeline output file of the pipeline
def generate_output_file(pipeline_data, pipeline_output_file, out_dir):

    """
    Generate the pipeline output file.

    Args:
        pipeline_data (PipelineData): Instance of the PipelineData class.
        pipeline_output_file (str): Name of the output file of the pipeline
        out_dir (str): Directory to save the output file.
    """

    if out_dir is None:
        out_dir = os.getcwd()
        print(out_dir)

    # Create the directory if it doesn't exist
    create_directory(out_dir) 

    # Get the contigs and bins from the PipelineData instance
    contigs = pipeline_data.get_contigs()
    bins = pipeline_data.get_bins()

    #the path to the file result of the pipeline
    output_file = os.path.join(out_dir, pipeline_output_file)

    #write in the file
    with open(output_file, "w") as f:
        f.write("Contigs:\n")
        for contig_id, sequence in contigs.items():
            f.write(f"{contig_id}: {sequence}\n")

        f.write("\nBins:\n")
        for bin_id, contig_names in bins.items():
            f.write(f"{bin_id}: {contig_names}\n")

    
    print(f"Results saved to {output_file}")

def main():


    try:
        args = parse_arguments()
        method_configs = load_yaml(args.user_file)

        # Get the pipeline input files
        input_gfa, input_fasta = get_input_file(method_configs)

        # Get pipeline configurations
        plaspipe_configs = get_pipeline_configs(method_configs)

        # Unpack the configurations
        classification_dir, binning_dir, prefix, class_tool_name, class_tool_version, class_input_format, bin_tool_name, bin_tool_version, bin_input_format, out_dir = plaspipe_configs

        print(f'Classification output directory: {classification_dir}')

        # Check if the pipeline can handle the input format
        if (class_input_format == 'GFA' or bin_input_format == 'GFA') and (input_gfa is None):
            print("This pipeline cannot handle the conversion of FASTA to GFA")
        else:
            # Create the PipelineData instance
            plaspipe_data = run_plaspipe_data(method_configs, prefix, class_tool_name, class_tool_version, bin_tool_name, bin_tool_version, input_gfa, input_fasta, classification_dir, binning_dir)

            # Generate the pipeline output file
            generate_output_file(plaspipe_data, args.pipeline_output, out_dir)

    except Exception as e:
        process_exception(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

