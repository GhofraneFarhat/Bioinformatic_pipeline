import yaml
from classification_wrapper import ClassificationWrapper
import os
from pipeline_data import PipelineData
from binning_wrapper import BinningWrapper
import json
import time
import argparse



# Create the argument parser
parser = argparse.ArgumentParser(description='Pipeline for classification and binning')

# Add arguments
parser.add_argument('user_file', type=str, help='path to the yaml file that the user with provide it to the pipeline')
parser.add_argument('pipeline_input_file', type=str, help='path to pipeline input file, this input file can be fasta or GFA or both ')
parser.add_argument('pipeline_output', type=str, help='the file that the contigs dict and bins dict will be saved ther, the output of the pipeline')

# Parse the command-line arguments
args = parser.parse_args()

# Access the arguments
user_file = args.user_file
pipeline_input_file = args.pipeline_input_file
pipeline_output = args.pipeline_output




# Load the YAML file
#user_file = 'user.yaml'
with open(user_file, 'r') as yaml_file:
    method_configs = yaml.safe_load(yaml_file)

#define the input path
#input_file = "output.fasta"
#pipeline_input_file = "binniginput.GFA"
input_gfa = ''
input_fasta = ''
if pipeline_input_file.endswith('.GFA'):
    input_gfa = pipeline_input_file
elif pipeline_input_file.endswith('.fasta'):
    input_fasta = pipeline_input_file


#get the work direct
current_dir = os.getcwd() #the director of work


#for classification
#get the classification keys
classification_config = method_configs['classification']#contains the name, the version, the command

#for binning
# save the binning configuration from the yaml file
binning_config = method_configs['binning'] #contains the name, the version, the command of the binning tool

#use time to generate output files
timestamp = int(time.time())  

#need to figure out about the conv.fasta
#conv_file = "conv.fasta"
conv_file = f"conv_{timestamp}.fasta" 

#create the pipline instance
pipeline_data = PipelineData(input_gfa, input_fasta, conv_file )

#condition for the format: if one of tools require a gfa as an input format and the pipeline input is fasta then we will get an error
if (classification_config['input_format'] == 'GFA' or binning_config['input_format'] == 'GFA') and (pipeline_input_file.endswith ('fasta')):
    print("this pipeline can't handle the conversion of fasta to GFA")

else:
    #load the contigs from the input file
    pipeline_data.load_contigs()




    #**********************************
    #classification part:

    # Define the output file
    #output_classification = "classification.fasta"
    output_classification = f"classification_{timestamp}.fasta"
    #there is an error when running if the inputclass(input file) isn't in the direct


    
    # get the path to classification tool
    repo_path = current_dir + '\\' + classification_config.get('name')

    #create the classwrapper instance
    classification_wrapper = ClassificationWrapper(output_classification, classification_config, repo_path, pipeline_data)

    #run the classwrapper 
    classification_wrapper.run()




    #********problems:

    #need to try the converted input_file
    #need to see about the output_file
    #*********************************


    ##############################################
    ##############################################

    #binning part


    #the output file for the binning
    #i want to generate it automatic

    #output_path = "result.csv"
    output_path = f"result_{timestamp}.csv"

    #get the binning tool repo path 
    bin_path = current_dir + '\\' + binning_config.get('name')

    #create the binwrapper instance
    binning_wrapper = BinningWrapper(input_gfa, output_path, binning_config, bin_path, pipeline_data)

    #run the binwrapper
    binning_wrapper.run()






    #************Problems:
    #the converted file 
    #if bin tool demands a GFA file and the pipeline input is a fasta just print an error
    #if the bin tool demands a fasta file and the pipeline input is a GFA file then convert


    #result of the pipeline 

    #get the result of the classification 
    contigs = pipeline_data.get_contigs()

    # get the contigs result from the binning tool
    bins = pipeline_data.get_bins()
    
    """
    for bin_id, contig_names in bins.items():
        
        print(f' {bin_id}: {contig_names}')

    
    
    
    for contig_id, sequence in contigs.items():
        print(f"{contig_id}: {sequence}")

    """
    #pipeline_output_file = "output.txt"
    output_file = os.path.join(current_dir, pipeline_output)
    with open(output_file, "w") as f:
            f.write("Contigs:\n")
            for contig_id, sequence in contigs.items():
                f.write(f"{contig_id}: {sequence}\n")

            f.write("\nBins:\n")
            for bin_id, contig_names in bins.items():
                f.write(f"{bin_id}: {contig_names}\n")

    print(f"Results saved to {output_file}")

