import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import sys
from .tool_command import get_command
from .tool_conversion import run_conversion
from .gfa_to_fasta import write_GFA_to_FASTA


class ClassificationWrapper:

    """
    the classification wrapper is reponsible for the formating the input of the pipeline (fasta or GFA) and for
    converting the output of the tool format to the pipeline output format (csv)

    Args: 
    Returns:
    """ 

#initialisation
    def __init__(self, classification_folder, prefix, method_config, fasta_path="", gfa_path = "", gzipped_gfa=False, gzipped_fasta=False):
        #self.input_path = input_file
        #self.output_classif = output_of_classification
        self.classification_dir = classification_folder
        self.prefix = prefix
        self.method_configuration = method_config
        #self.repo_path = repo_path
        #self.output_file = output_conv
        self.fasta_path = fasta_path
        self.gfa_path = gfa_path

        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta
        #self.id_fun = id_fun
        #self.class_dir = temp_dir
        #self.pipeline_data = pipeline_data
        #the path to save the output of the tool 
        #self.class_conversion = class_conv  






    #to run the tool 
    def run(self):


        # Run the classification tool : need to see more about the parameters
        resultat_of_wrapper = self.run_classification()
        
        return resultat_of_wrapper


    #the function to run the classification tool
    def run_classification(self):
        
        """ 
        args:
        command line: usually python script_classification.py
        parameters: the tool parameter
        input_classif_converted: the input file for the classifcation tool
        output_classif: the output file result of running the classification tool

        return: 
        resultat the classification tool
        """ 


        #extract the input format
        class_format = self.method_configuration['input_format']  # class_format for the format of the input of the classification tool
        
        #extract informations about the classification tool
        class_tool_name = self.method_configuration['name']
        #prefix = self.method_configuration['prefix']
        version = self.method_configuration['version']

        #extarct tool format of the output: need to find if i can do it without
        output_format = self.method_configuration['output_format']

        # Absolute path to the tool directory
        tool_directory = 'C:/Users/user/Desktop/Bioinformatic_pipeline/plASgraph'

        # Change the working directory to the tool directory
        os.chdir(tool_directory)

        #Formatting the input file of the pipeline 
        self.input_classif_converted = self.conversion(class_format, class_tool_name) #input_classif_converted: the path to input file
        print(f'here is my converted file {self.input_classif_converted}')
        #input_classif_converted = 'C:/Users/user/Desktop/run/input/output.fasta'

        # getting the tool parameter using the yamk file
        tool_parameters = self.method_configuration['parameters'] # the parameter of the classification tool

        #the file that the output results of the classification tool will be saved, have the format of the classification_tool
        output_classification = os.path.join(self.classification_dir, self.prefix + '_' + class_tool_name + '_' + version + '.' + output_format )
        


        # Construct the full command with input_file_path and output_file
        #full_command = command_parts + ['-i', self.input_classif_converted, '-o', output_classification]
        full_command = get_command(self.method_configuration, self.input_classif_converted, output_classification)
        # Run the command using subprocess
        subprocess.run(full_command, shell=False)

        #need to return the output_classif
        return output_classification


    #convert to csv file
    def pipeline_conversion_to_csv(self, file_path):
        
        """ 
        this is the place of the conversion: 
        args: 
        tool_name
        tool_version
        path to tool_result_file
        return:
        pipeline_file: the standart file that we will use for the pipeline to update contigs 

        """ 

        #extract informations about the classification tool
        class_tool_name = self.method_configuration['name']
        #prefix = self.method_configuration['prefix']
        version = self.method_configuration['version']
        


        #specify csv file path
        csv_file = csv_file = os.path.join(self.classification_dir, self.prefix + '_' + class_tool_name + '_' + version + '.csv')
        print(f'this is my csv file enjoy {csv_file}')


        # Run the conversion script with run conversion
        run_conversion(class_tool_name, version, file_path, csv_file)
        print(f'i am the csv file that u are searching for me {csv_file}')
        return csv_file

    #getter for the classification_wrapper
    def get_csv_file(self):
        
        """ 
        returns: the csv file for the pipeline
        """
        resultat_of_classification = self.run()

        output_pipeline_csv = self.pipeline_conversion_to_csv(resultat_of_classification)
        print("hey now I'm a csv file")
        return output_pipeline_csv 


    ######################################################################
    ######################################################################
    #conversion from GFA to fasta
    

    # The conversion function that will provide the path to the wrapper
    #the conversionn function
    def conversion(self, input_format, folder_tool):

        """
        conversion of gfa to fasta or just return path

        Args:
        input_format: fasta or gfa
        folder_tool: to create the output file on the desired directory
        gfa_path: the path to the gfa file, input of the pipeline
        fasta_path: the path to the fasta file, input of the pipeline
        gzipped: file gzipped or not
        output_file: where the converted input got saved

        return:
            the path to the same file path if there's no conversion
            the path to the converted file 

        """
        #print(self.output_file)
        print("helllllllllllllllllllllllo")
        print(self.gfa_path)

        # If classification_dir is None, set it to the current working directory
        if self.classification_dir is None:
            self.classification_dir = os.getcwd()
        
        #prefix = self.method_configuration['prefix']
        output_classif_file = os.path.join (self.classification_dir , self.prefix + '_converted.fasta')
        print(f"hello i wanna be converted here {output_classif_file}")

        # Create the output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_classif_file), exist_ok=True)
        
        #classify_dir_output = os.path.join(os.getcwd() + '\\' + folder_tool , self.output_file)

        #print(f'this the dir{classify_dir_output}')

        sep = ' '
        if input_format == 'fasta':
            if self.fasta_path:
                return self.fasta_path
            elif self.gfa_path:
                print("wanna convert gfa with me?")
                #problem with the the tool name 
                #self.write_GFA_to_FASTA(self.gfa_path, working_dir_output, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                write_GFA_to_FASTA(self.gfa_path, output_classif_file, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                print(f"let's convert this gfa file {output_classif_file}")
                return output_classif_file
                
            else:
                raise ValueError("Neither GFA nor FASTA file provided as input.")
        elif input_format == 'gfa' or input_format == 'GFA':
            if self.gfa_path:
                return self.gfa_path
            else:
                raise ValueError("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
        else:
            raise ValueError(f"You have to provide a fasta or GFA, Unsupported input format: {input_format}")
        

