import yaml
from classification_wrapper import ClassificationWrapper
import os
from pipeline_data import PipelineData

# Load the YAML file
with open('user.yaml', 'r') as yaml_file:
    method_configs = yaml.safe_load(yaml_file)

#define the input path
input_file = "input_pipeline.fasta"

#create the pipline instance
pipeline_data = PipelineData(input_file)

#load the contigs from the input file
pipeline_data.load_contigs()

#get the work direct
current_dir = os.getcwd() #the director of work

#get the classification keys
classification_config = method_configs['classification']#contains the name, the version, the command
 


#**********************************
#classification part:
# Define the input file path and output file
output_file = "res.fasta"
#there is an error when running if the inputclass(input file) isn't in the direct

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