import yaml
from .classification_wrapper import ClassificationWrapper
import os
from .plaspipe_data import PipelineData
from .binning_wrapper import BinningWrapper
import json
import time
import argparse



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
    gfa_path = method_config['input']['path_to_input_gfa']
    fasta_path = method_config['input']['path_to_input_fasta']
    return gfa_path, fasta_path



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

    #create the instance of plaspipe_data
    plaspipe_data = PipelineData(classification_wrapper, class_name, class_version, input_gfa, input_fasta)

    #load contigs name
    plaspipe_data.load_contigs()

    #run the classification wrapper 
    classification_wrapper.run()

    #update the plaspipe_data
    plaspipe_data.update_plaspipe_data_from_classwrapper()

    #return the plaspipe_data after the changement of the classification wrapper, add contigs class
    return plaspipe_data

#the function to run the binning wrapper
def run_binning(config, plaspipe_data, input_gfa, timestamp):
    output_path = f"result_{timestamp}.csv"
    binning_wrapper = BinningWrapper(input_gfa, output_path, config, plaspipe_data)
    binning_wrapper.run()



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
    prefix = method_configs['prefix']
    #get the work direct
    current_dir = os.getcwd()
    

    #use time to generate output files
    timestamp = int(time.time())



    if (method_configs['classification']['input_format'] == 'GFA' or method_configs['binning']['input_format'] == 'GFA') and (input_gfa == None):
        print("this pipeline can't handle the conversion of fasta to GFA")
    else:

        #the return of the run_classification function is pipeline data
        plaspipe_data = run_classification(classification_dir, prefix, method_configs['classification'], input_fasta, input_gfa,method_configs['classification']['name'],method_configs['classification']['version'] )

        #binning_repo_path = os.path.join(current_dir, method_configs['binning']['name'])
        #run_binning(method_configs['binning'], plaspipe_data, input_gfa, timestamp)

        #for the output of the pipeline
        out_dir = method_configs['outdir_pipeline']
        generate_output_file(plaspipe_data, args.pipeline_output, out_dir)

if __name__ == "__main__":
    main()

