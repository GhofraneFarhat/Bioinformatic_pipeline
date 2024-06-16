import yaml
from .classification_wrapper import ClassificationWrapper
import os
from .plaspipe_data import PipelineData
from .binning_wrapper import BinningWrapper
import json
import time
import argparse
import gzip
import sys
import shutil
import logging



def parse_arguments():
    """ 
    reading the argument from the command line
    return: args: list of arguments
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(description='Pipeline for bioinformatic tools (classification and binning tools for plasmid analysis')

    # Add arguments
    parser.add_argument('user_file', type=str, help='path to the yaml file that the user with provide it to the pipeline')
    #parser.add_argument('pipeline_input_file', type=str, help='path to pipeline input file, this input file can be fasta or GFA or both ')
    parser.add_argument('pipeline_output', type=str, help='the file that the contigs dict and bins dict will be saved ther, the output of the pipeline')

    # Parse the command-line arguments
    args = parser.parse_args()

    return args


# the yaml file, from the user need some changes
def load_yaml(file):
    """ load the yaml file """
    """ 
    Args: file: the yaml file
    Returns: contigs
    """
    with open(file, 'r') as yaml_file:
        configs = yaml.safe_load(yaml_file)
    return configs


#the input file of the pipeline can be fasta or gfa or both
#need to change, the inpute of the pipeline will be provided from the yaml file
def get_input_file(method_config):

    outdir_pipeline = method_config['outdir_pipeline']
    os.makedirs(outdir_pipeline, exist_ok=True)

    gun_gfa_path = os.path.join(method_config['outdir_pipeline'], "gunzipped_gfa.gfa") 
    gun_fasta_path = os.path.join(method_config['outdir_pipeline'], "gunzipped_fasta.fasta")  

    gfa_path = method_config['input']['path_to_input_gfa']
    fasta_path = method_config['input']['path_to_input_fasta']
    print(f"this is my gfa file {gfa_path}")

    if gfa_path and gfa_path.endswith('.gz'):
        gfa_path = gunzip_GFA(gfa_path, gun_gfa_path)
        
            
    if fasta_path and fasta_path.endswith('.gz'):
        fasta_path = gunzip_FASTA(fasta_path, gun_fasta_path)

    return gfa_path, fasta_path

#gunzzipping a gfa file
def gunzip_GFA(in_file_path, out_file_path):
    """
    Gunzip a GFA file

    Args:
       in_file_path (str): path to input gzipped GFA file
       out_file_path (str): path to output GFA file

    Returns:
       Creates GFA file out_file_path
    """
    try:
        with gzip.open(in_file_path) as in_file, open(out_file_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
            print(f'this my gfa file gunzipped {out_file_path}')
            return out_file_path
    except Exception as e:
        process_exception(f'FASTA\tGunzipping {in_file_path} to {out_file_path}: {e}')

#gunzip a fasta file
def gunzip_FASTA(in_file_path, out_file_path):
    """
    Gunzip a FASTA file

    Args:
       in_file_path (str): path to input gzipped FASTA file
       out_file_path (str): path to output FASTA file

    Returns:
       Creates FASTA file out_file_path
    """
    records = []
    with gzip.open(in_file_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
        with open(out_file_path, 'w') as out_file:
            SeqIO.write(records, out_file, 'fasta')

def process_exception(msg):
    logging.exception(msg)
    print(f'EXCEPTION\t{msg}', file=sys.stderr)
    sys.exit(1)




#the function to run the classification wrapper
def run_classification(classification_folder, prefix, config, input_fasta, input_gfa, class_name, class_version):

    """ run the classifcation wrapper"""
    """ 
    Args: 
    config: the configuration of the classification tool
    repo_path : the path of the classification tool folder
    conv_file: the conversion file, the result of the conversion function
    input_fasta, input_gfa: the input of the pipeline
    Returns:
    plaspipe_data: the data of the pipeline updated after passing by the classification tool
    """


    #create the classwrapper instance
    classification_wrapper = ClassificationWrapper(classification_folder, prefix, config, input_fasta, input_gfa)

    #run the classification wrapper 
    classification_wrapper.run()


    #return the plaspipe_data after the changement of the classification wrapper, add contigs class
    return classification_wrapper

#the function to run the binning wrapper
def run_binning(binning_folder, prefix, config, input_fasta, input_gfa, classification_path, class_tool_name, class_tool_version):
    
    print(f'classification folder {classification_path}')
    binning_wrapper = BinningWrapper(binning_folder, prefix, config, input_fasta, input_gfa, classification_path, class_tool_name, class_tool_version)

    binning_wrapper.run()

    return binning_wrapper

def run_plaspipe_data (classification_wrapper, binning_wrapper, class_tool_name, class_tool_version, input_gfa, input_fasta):

    #cretae the plaspipe_data instance
    plaspipe_data = PipelineData(classification_wrapper, binning_wrapper, class_tool_name, class_tool_version, input_gfa, input_fasta)

    #load contigs name
    plaspipe_data.load_contigs()

    #update the plaspipe_data
    plaspipe_data.update_plaspipe_data_from_classwrapper()
        
    #update the plaspipe_data
    plaspipe_data.update_plaspipe_data_from_binwrapper()

    return plaspipe_data




#the function to generate the pipeline output file of the pipeline
def generate_output_file(pipeline_data, pipeline_output_file, out_dir):
    """ 
    this function is to generate the pipeline result saved in a file
    Args:
    pipeline_data: the flow of data that is updated all over the pipeline
    pipeline_output_file: the output file for the command line"bin_tool
    temp_dir: the path where the output file will be saved
    """
    if out_dir is None:
        out_dir = os.getcwd()
        print(out_dir)

    # Create the directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)

    #the result of the classification tool
    contigs = pipeline_data.get_contigs()

    #the result of the binning tool
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
    args = parse_arguments()

    method_configs = load_yaml(args.user_file)

    #the pipeline input file
    #need to figure out if the user provide two inputs (fa and gfa)
    input_gfa, input_fasta = get_input_file(method_configs)

    #extract the temp directory of the output different format of the classification tool
    classification_dir = method_configs['classification']['output']['outdir_classification']
    binning_dir = method_configs['binning']['output']['outdir_binning']
    prefix = method_configs['prefix']
    
    print(f'the classifcation_dir {classification_dir}')
    
    
    
    class_tool_name = method_configs['classification']['name']
    class_tool_version = method_configs['classification']['version']




    if (method_configs['classification']['input_format'] == 'GFA' or method_configs['binning']['input_format'] == 'GFA') and (input_gfa == None):
        print("this pipeline can't handle the conversion of fasta to GFA")
    else:

        #create the instance of plaspipe_data
        #plaspipe_data = PipelineData(class_tool_name, class_tool_version, input_gfa, input_fasta)
        



        #run the classiifcation wrapper function
        classification_wrapper = run_classification(classification_dir, prefix, method_configs['classification'], input_fasta, input_gfa, class_tool_name, class_tool_version )

        #run the binning wrapper function
        binning_wrapper = run_binning(binning_dir, prefix, method_configs['binning'], input_fasta, input_gfa, classification_dir, class_tool_name, class_tool_version)

        #run the plaspipe function
        plaspipe_data = run_plaspipe_data(classification_wrapper, binning_wrapper, class_tool_name, class_tool_version, input_gfa, input_fasta)

        #for the output of the pipeline
        out_dir = method_configs['outdir_pipeline']
        generate_output_file(plaspipe_data, args.pipeline_output, out_dir)

if __name__ == "__main__":
    main()

