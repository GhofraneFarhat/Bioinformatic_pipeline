import yaml
from classification_wrapper import ClassificationWrapper
import os
from plaspipe_data import PipelineData
from binning_wrapper import BinningWrapper
import json
import time
import argparse

def read_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Pipeline for classification and binning')

    # Add arguments
    parser.add_argument('user_file', type=str, help='path to the yaml file that the user with provide it to the pipeline')
    parser.add_argument('pipeline_input_file', type=str, help='path to pipeline input file, this input file can be fasta or GFA or both ')
    parser.add_argument('pipeline_output', type=str, help='the file that the contigs dict and bins dict will be saved ther, the output of the pipeline')

    # Parse the command-line arguments
    args = parser.parse_args()

    return args

def read_yaml_file(user_file):
    # Load the YAML file
    with open(user_file, 'r') as yaml_file:
        method_configs = yaml.safe_load(yaml_file)

    return method_configs

def run_classification_wrapper(classification_config, repo_path, pipeline_data, current_dir, timestamp):
    # Define the output file
    output_classification = f"classification_{timestamp}.fasta"

    # Create the ClassificationWrapper instance
    classification_wrapper = ClassificationWrapper(output_classification, classification_config, repo_path, pipeline_data)

    # Run the classification wrapper
    classification_wrapper.run()

def run_binning_wrapper(binning_config, bin_path, pipeline_data, input_gfa, timestamp):
    # Define the output file
    output_path = f"result_{timestamp}.csv"

    # Create the BinningWrapper instance
    binning_wrapper = BinningWrapper(input_gfa, output_path, binning_config, bin_path, pipeline_data)

    # Run the binning wrapper
    binning_wrapper.run()

def generate_output_file(pipeline_data, pipeline_output, current_dir):
    # Get the result of the classification
    contigs = pipeline_data.get_contigs()

    # Get the contigs result from the binning tool
    bins = pipeline_data.get_bins()

    output_file = os.path.join(current_dir, pipeline_output)
    with open(output_file, "w") as f:
        f.write("Contigs:\n")
        for contig_id, sequence in contigs.items():
            f.write(f"{contig_id}: {sequence}\n")

        f.write("\nBins:\n")
        for bin_id, contig_names in bins.items():
            f.write(f"{bin_id}: {contig_names}\n")

    print(f"Results saved to {output_file}")

def main():
    # Read command-line arguments
    args = read_arguments()

    # Access the arguments
    user_file = args.user_file
    pipeline_input_file = args.pipeline_input_file
    pipeline_output = args.pipeline_output

    # Get the current working directory
    current_dir = os.getcwd()

    # Read the YAML file
    method_configs = read_yaml_file(user_file)

    # Get the input file formats
    input_gfa = ''
    input_fasta = ''
    if pipeline_input_file.endswith('.GFA'):
        input_gfa = pipeline_input_file
    elif pipeline_input_file.endswith('.fasta'):
        input_fasta = pipeline_input_file

    # Get the classification configuration
    classification_config = method_configs['classification']

    # Get the binning configuration
    binning_config = method_configs['binning']

    # Generate a timestamp for output files
    timestamp = int(time.time())

    # Create the conversion file name
    conv_file = f"conv_{timestamp}.fasta"

    # Create the PipelineData instance
    pipeline_data = PipelineData(input_gfa, input_fasta, conv_file)

    # Check if the input format is compatible with the tools
    if (classification_config['input_format'] == 'GFA' or binning_config['input_format'] == 'GFA') and pipeline_input_file.endswith('fasta'):
        print("This pipeline can't handle the conversion of fasta to GFA")
    else:
        # Load the contigs from the input file
        pipeline_data.load_contigs()

        # Get the path to the classification tool
        repo_path = os.path.join(current_dir, classification_config.get('name'))

        # Run the classification wrapper
        run_classification_wrapper(classification_config, repo_path, pipeline_data, current_dir, timestamp)

        # Get the path to the binning tool
        bin_path = os.path.join(current_dir, binning_config.get('name'))

        # Run the binning wrapper
        run_binning_wrapper(binning_config, bin_path, pipeline_data, input_gfa, timestamp)

        # Generate the output file
        generate_output_file(pipeline_data, pipeline_output, current_dir)

if __name__ == "__main__":
    main()