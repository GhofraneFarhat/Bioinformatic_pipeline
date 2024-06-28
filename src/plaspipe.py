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
from .plaspipe_utils import process_error
from .plaspipe_utils import check_file
from .plaspipe_utils import gunzip_FASTA
from .plaspipe_utils import gunzip_GFA
from .plaspipe_utils import create_directory
from .plaspipe_utils import verify_input_file
from .plaspipe_utils import check_gfa_input
from .plaspipe_utils import process_file
from .plaspipe_utils import verif_file
from .plaspipe_utils import exist_file
from .plaspipe_utils import setup_logging 
from .plaspipe_utils import log_file_creation


from .plaspipe_data import PipelineData



import os
import argparse

import os
import argparse

def parse_arguments():
    """
    Read arguments from the command line.
    Returns:
        args (Namespace): Parsed arguments.
    """

    try:
        # Create the argument parser
        parser = argparse.ArgumentParser(description='Pipeline for bioinformatic tools (classification and binning tools for plasmid analysis')

        # Add optional arguments with flags
        parser.add_argument('-yf', '--user_file', type=str, required=True, help='Path to the YAML file provided by the user')
        #parser.add_argument('-out', '--pipeline_output', type=str, required=True, help='File to save contigs dict and bins dict (pipeline output)')

        # Parse the arguments
        args = parser.parse_args()

        # Verify the user_file is a YAML file
        verif_file(args.user_file, '.yaml')

        # Check if the YAML file exists
        exist_file(args.user_file)

        # Check if the output directory exists
        #exist_file(args.pipeline_output)

        return args

    except argparse.ArgumentError as ae:
        process_error(f"Argument parsing error: {ae}")
        sys.exit(1)
    except argparse.ArgumentTypeError as ate:
        process_error(f"Argument type error: {ate}")
        sys.exit(1)




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
            raise ValueError(f'YAML file {file} is empty.')
        return configs
    except FileNotFoundError as f:
        process_file(f"The file is not found: {f}")





def get_input_file(method_config):

    outdir_pipeline = method_config['outdir_pipeline']
    prefix = method_config['prefix']

    pipeline_output_path = create_directory(outdir_pipeline, prefix)

    gun_gfa_path = os.path.join(pipeline_output_path, "gunzipped_gfa.gfa") 
    gun_fasta_path = os.path.join(pipeline_output_path, "gunzipped_fasta.fasta")  

    gfa_path = method_config['input']['path_to_input_gfa']
    fasta_path = method_config['input']['path_to_input_fasta']
    print(f"this is my gfa file {gfa_path}")
    log_file_creation('gfa_path', gfa_path)
    log_file_creation('fasta_path', fasta_path)

    #unzip the gfa file
    try:
        logging.info(f"Processing GFA file: {gfa_path}")
        if gfa_path and gfa_path.endswith('.gz'):
            gfa_path = gunzip_GFA(gfa_path, gun_gfa_path)
    except Exception as e:
        logging.error(f"Error unzipping GFA file {gfa_path}: {e}")


    #unzip the fasta file
    try:
        logging.info(f"Processing FASTA file: {fasta_path}")
        if fasta_path and fasta_path.endswith('.gz'):
            fasta_path = gunzip_FASTA(fasta_path, gun_fasta_path)
    except Exception as e:
        logging.error(f"Error unzipping FASTA file {fasta_path}: {e}")

    #chech the input files
    verify_input_file(gfa_path, fasta_path)


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
    pipeline_configs.extend([class_tool_name, class_tool_version, class_input_format, bin_tool_name, bin_tool_version, bin_input_format])

    # Extract output directory for the pipeline
    out_dir = method_configs['outdir_pipeline']
    pipeline_configs.append(out_dir)

    return pipeline_configs


#the function to generate the pipeline output file of the pipeline
import csv

def generate_output_files(pipeline_data, out_dir, prefix):

    """
    Generate the pipeline output files

    Args:
        pipeline_data (PipelineData): Instance of the PipelineData class
        out_dir (str): Directory to save the output files
    """

    # Create the directory if it doesn't exist
    out_dir = create_directory(out_dir, prefix) 

    # Get the contigs and bins from the PipelineData instance
    contigs = pipeline_data.get_contigs()
    bins = pipeline_data.get_bins()

    # Define the paths to the output files
    binning_output_file = os.path.join(out_dir, "binning_tool_output.csv")
    classification_output_file = os.path.join(out_dir, "classification_tool_output.csv")

    # Write the contigs to a CSV file
    with open(classification_output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Contig ID", "Sequence"])
        for contig_id, sequence in contigs.items():
            writer.writerow([contig_id, sequence])

    log_file_creation("classification_output_file", classification_output_file)

    # Write the bins to a CSV file
    with open(binning_output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Bin ID", "Contig Names"])
        for bin_id, contig_names in bins.items():
            writer.writerow([bin_id, ', '.join(contig_names)])

    log_file_creation("binning_output_file", binning_output_file)

    print(f"Binning results saved to {binning_output_file}")
    print(f"Classification results saved to {classification_output_file}")


def main():
    try:
        args = parse_arguments()
        method_configs = load_yaml(args.user_file)
        setup_logging(method_configs['outdir_pipeline'])
        # Get the pipeline input files
        input_gfa, input_fasta = get_input_file(method_configs)

        # Get pipeline configurations
        plaspipe_configs = get_pipeline_configs(method_configs)

        # Unpack the configurations
        classification_dir, binning_dir, prefix, class_tool_name, class_tool_version, class_input_format, bin_tool_name, bin_tool_version, bin_input_format, out_dir = plaspipe_configs

        print(f'Classification output directory: {classification_dir}')

        # Check if the pipeline can handle the input format
        try:

            #verify if the input format of the tools is supported
            check_gfa_input(class_input_format, bin_input_format, input_gfa)

            # Create the PipelineData instance
            plaspipe_data = run_plaspipe_data(method_configs, prefix, class_tool_name, class_tool_version, bin_tool_name, bin_tool_version, input_gfa, input_fasta, classification_dir, binning_dir)

            # Generate the pipeline output file
            generate_output_files(plaspipe_data, out_dir, prefix)

        except ValueError as ve:
            process_error(str(ve))

    except Exception as e:
        process_exception(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

