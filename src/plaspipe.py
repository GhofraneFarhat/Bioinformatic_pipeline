import yaml
from .classification_wrapper import ClassificationWrapper
import os
from .plaspipe_data import PipelineData
from .binning_wrapper import BinningWrapper
import json
import time
import argparse


#reading the arguments from the command line
def parse_arguments():
    # Create the argument parser
    parser = argparse.ArgumentParser(description='Pipeline for bioinformatic tools (classification and binning tools for plasmid analysis')

    # Add arguments
    parser.add_argument('user_file', type=str, help='path to the yaml file that the user with provide it to the pipeline')
    parser.add_argument('pipeline_input_file', type=str, help='path to pipeline input file, this input file can be fasta or GFA or both ')
    parser.add_argument('pipeline_output', type=str, help='the file that the contigs dict and bins dict will be saved ther, the output of the pipeline')

    # Parse the command-line arguments
    args = parser.parse_args()

    return args


#load the yaml file, from the user and need some changes
def load_yaml(file):
    with open(file, 'r') as yaml_file:
        configs = yaml.safe_load(yaml_file)
    return configs


#the input file of the pipeline can be fasta or gfa or both
def get_input_file_type(file):
    if file.endswith('.GFA'):
        return file, ''
    elif file.endswith('.fasta'):
        return '', file


#the function to run the classification wrapper
def run_classification(config, repo_path, conv_file, input_fasta, input_gfa, timestamp):

    """ run the classifcation wrapper"""
    """ 
    Args: 
    config: the configuration of the classification tool
    repo_path : the path of the classification tool folder
    conv_file: the conversion file, the result of the conversion function
    """
    #the file that the output results of the classification tool will be saved, have the format of the classification_tool
    output_classification = f"classification_{timestamp}.fasta"
    
    #create the classwrapper instance
    classification_wrapper = ClassificationWrapper(output_classification, config, repo_path, conv_file, input_fasta, input_gfa)

    #create the instance of plaspipe_data
    pipeline_data = PipelineData(classification_wrapper, input_gfa, input_fasta)

    #load contigs name
    pipeline_data.load_contigs()

    #run the classification wrapper 
    classification_wrapper.run()

    #update the plaspipe_data
    pipeline_data.update_pipeline_data_from_classwrapper()

    #return the plaspipe_data after the changement of the classification wrapper, add contigs class
    return pipeline_data

def run_binning(config, bin_path, pipeline_data, input_gfa, timestamp):
    output_path = f"result_{timestamp}.csv"
    binning_wrapper = BinningWrapper(input_gfa, output_path, config, bin_path, pipeline_data)
    binning_wrapper.run()

def generate_output_file(pipeline_data, pipeline_output, current_dir):

    contigs = pipeline_data.get_contigs()

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
    args = parse_arguments()

    method_configs = load_yaml(args.user_file)

    #the pipeline input file
    #need to figure out if the user provide two inputs (fa and gfa)
    input_gfa, input_fasta = get_input_file_type(args.pipeline_input_file)

    #get the work direct
    current_dir = os.getcwd()

    #use time to generate output files
    timestamp = int(time.time())

    #the conversion files of the conversion of gfa to fasta
    conv_file = f"conv_{timestamp}.fasta"

    if (method_configs['classification']['input_format'] == 'GFA' or method_configs['binning']['input_format'] == 'GFA') and (args.pipeline_input_file.endswith ('fasta')):
        print("this pipeline can't handle the conversion of fasta to GFA")
    else:

        # get the path to classification tool
        classification_repo_path = os.path.join(current_dir, method_configs['classification']['name'])

        pipeline_data = run_classification(method_configs['classification'], classification_repo_path, conv_file, input_fasta, input_gfa, timestamp)

        binning_repo_path = os.path.join(current_dir, method_configs['binning']['name'])
        run_binning(method_configs['binning'], binning_repo_path, pipeline_data, input_gfa, timestamp)

        generate_output_file(pipeline_data, args.pipeline_output, current_dir)

if __name__ == "__main__":
    main()