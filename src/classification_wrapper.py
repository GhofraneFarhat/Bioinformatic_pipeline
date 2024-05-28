import os
import subprocess
from Bio import SeqIO


class ClassificationWrapper:


#initialisation
    def __init__(self, output_file, method_config, repo_path, pipeline_data):
        #self.input_path = input_file
        self.output_classif = output_file
        self.method_configuration = method_config
        self.repo_path = repo_path
        #self.class_dir = temp_dir
        self.pipeline_data = pipeline_data
        #the path to save the output of the tool        


#to run the tool frol the file load.py
    def run(self):


        # Run the classification tool : need to see more about the parameters
        resultat_of_classification = self.run_classification()

        # Update the PipelineData instance with the results
        self.update_pipeline_data(resultat_of_classification)


#the function to run the classification tool
    def run_classification(self):

        """ 
        args:
        command line: usually python script_classification.py
        parameters: the tool parameter
        input_classif: the input file for the classifcation tool
        output_classif: the output file result of running the classification tool

        return: 
        resultat the classification tool
        """ 
        # Extract the tool_command from the method_config dictionary
        tool_command = self.method_configuration.get("command")

        #extract the input format
        class_format = self.method_configuration['input_format']  # class_format for the format of the input of the classification tool
        
        #extract the directory of the classification tool
        class_tool = self.method_configuration['name']

        #Formatting the input file of the pipeline 
        self.input_classif = self.pipeline_data.conversion(class_format, class_tool) #input_classif: the path to input file
        
        # getting the tool parameter using the yamk file
        tool_parameters = self.method_configuration['parameters'] # the parameter of the classification tool

        # Split the tool_command into a list
        command_parts = tool_command.split()

        # Construct the full command with input_file_path and output_file
        full_command = command_parts + [(self.input_classif), self.output_classif]

        # Change the current working directory to the repository path
        os.chdir(self.repo_path)

        # Run the command using subprocess
        subprocess.run(full_command, shell=False)

        #need to return the output_classif
        return self.output_classif


#the function to update the pipeline instance

    # update the pipeline data using the classification tool
    def update_pipeline_data(self, result_of_classification_tool):
        print("let's update")
        #output_class = os.path.join(self.class_dir, 'class.fasta')
        output_class = result_of_classification_tool
        #output_class = 'classification.fasta'

        with open(output_class, 'r') as classification_result:
            print(f"I'm using {output_class}")
            #output csv
            if self.method_configuration['output_format'] == 'csv':
                self.update_pipeline_data_from_csv(classification_result)
            #output tsv
            elif self.method_configuration['output_format'] == 'tsv':
                self.update_pipeline_data_from_tsv(classification_result)
            #output txt
            elif self.method_configuration['output_format'] == 'txt':
                self.update_pipeline_data_from_txt(classification_result)
            #output fasta
            elif self.method_configuration['output_format'] == 'fasta':
                self.update_pipeline_data_from_fasta(classification_result)
            else:
                print("Unsupported output format")

                #we can add other tools or find a solution for a general format of the output

# update from a csv file
    def update_pipeline_data_from_csv(self, classification_result):
        for line in classification_result:
            
            contig_name, plasmid_score, chromosome_score = line.strip().split(',')
            contig_class = {
                'plasmid_score': float(plasmid_score),
                'chromosome_score': float(chromosome_score)
            }#dict of contigs from the classification tool
            self.pipeline_data.set_contigs(contig_name, contig_class)

    #we need to handle if in the conversation of the file we lose a contig cause of an error 


 # update from a tsv file 
    def update_pipeline_data_from_tsv(self, classification_result):
        for line in classification_result:
            contig_name, plasmid_score, chromosome_score = line.strip().split('\t')
            contig_class = {
                'plasmid_score': float(plasmid_score),
                'chromosome_score': float(chromosome_score)
            }
            self.pipeline_data.set_contigs(contig_name, contig_class)

# update from a txt file
    def update_pipeline_data_from_txt(self, classification_result):
        
        pass
#update from fasta file
    def update_pipeline_data_from_fasta(self, classification_result):
        
        for seq_record in SeqIO.parse(classification_result, "fasta"):
            contig_name = seq_record.id
            description = seq_record.description.split()
            #print(description[1])
            #need to find how to update from fasta
            if description[1] == 'chromosome':
                plasmid_score = 0
                chromosome_score = 1
                
            elif description[1] == 'plasmid':
                plasmid_score = 1
                chromosome_score = 0
            elif description[1] == '(ambiguous)':
                plasmid_score = 1
                chromosome_score = 1
                
            else:
                plasmid_score = 0
                chromosome_score = 0

            contig_class = {
                'plasmid_score': plasmid_score,
                'chromosome_score': chromosome_score
            }
            self.pipeline_data.set_contigs(contig_name, contig_class)
        #we need to handle if in the conversation of the file we lose a contig cause of an error 