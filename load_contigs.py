import yaml
from classification_wrapper import ClassificationWrapper
import os
from pipeline_data import PipelineData
from binning_wrapper import BinningWrapper



# Load the YAML file
with open('user.yaml', 'r') as yaml_file:
    method_configs = yaml.safe_load(yaml_file)

#define the input path
input_file = "input_pipeline.fasta"
input_path = "binniginput.GFA"

#create the pipline instance
pipeline_data = PipelineData(input_file)

#load the contigs from the input file
pipeline_data.load_contigs()

#get the work direct
current_dir = os.getcwd() #the director of work



#**********************************
#classification part:
# Define the output file
output_file = "res.fasta"
#there is an error when running if the inputclass(input file) isn't in the direct

#get the classification keys
classification_config = method_configs['classification']#contains the name, the version, the command
 
# get the path to classification tool
repo_path = current_dir + '\\' + classification_config.get('name')

#create the classwrapper instance
classification_wrapper = ClassificationWrapper(input_file, output_file, classification_config, repo_path, pipeline_data)

#run the classwrapper 
classification_wrapper.run()

#print the result of the classification 
contigs = pipeline_data.get_contigs()
print("Loaded contigs:")
for contig_id, sequence in contigs.items():
    print(f"{contig_id}: {sequence}")


#********problems:

#need to try the converted input_file
#need to see about the output_file
#*********************************


##############################################
##############################################

#binning part

# save the binning configuration from the yaml file
binning_config = method_configs['binning'] #contains the name, the version, the command of the binning tool

#the output file for the binning
output_path = "result.csv"

#get the binning tool repo path 
bin_path = current_dir + '\\' + binning_config.get('name')

#create the binwrapper instance
binning_wrapper = BinningWrapper(input_path, output_path, binning_config, bin_path, pipeline_data)

#run the binwrapper
binning_wrapper.run()

# print the contigs result from the binning tool
bins = pipeline_data.get_bins()
for bin_id, contig_names in bins.items():
    #print(f'Bin {bin_id}: {", ".join(contig_names)}')
    print(f' {bin_id}: {contig_names}')




#************Problems:
#the converted file 
#if bin tool demands a GFA file and the pipeline input is a fasta just print an error
#if the bin tool demands a fasta file and the pipeline input is a GFA file then convert




