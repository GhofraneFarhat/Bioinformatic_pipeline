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


class BinningWrapper:

    """
    the binning wrapper is reponsible for the formating the input of the pipeline (fasta or GFA) and for
    converting the output of the tool format to the pipeline output format (csv)

    Args: 
    Returns:
    """ 

#initialisation
    def __init__(self, binning_folder, prefix, binning_config, fasta_path="", gfa_path = "", gzipped_gfa=False, gzipped_fasta=False, id_fun=lambda x: x):
        #self.input_path = input_file
        #self.output_classif = output_of_classification
        self.binning_dir = binning_folder
        self.prefix = prefix
        self.method_configuration = binning_config
        #self.repo_path = repo_path
        #self.output_file = output_conv
        self.fasta_path = fasta_path
        self.gfa_path = gfa_path

        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta
        self.id_fun = id_fun
        #self.class_dir = temp_dir
        #self.pipeline_data = pipeline_data
        #the path to save the output of the tool 
        #self.class_conversion = class_conv  

        # Mandatory fields in GFA contigs and links
        self.GFA_SEQ_KEY = 'Sequence'
        self.GFA_LEN_KEY = 'Length'
        self.GFA_FROM_KEY = 'From'
        self.GFA_FROM_ORIENT_KEY = 'FromOrient'
        self.GFA_TO_KEY = 'To'
        self.GFA_TO_ORIENT_KEY = 'ToOrient'
        self.GFA_OVERLAP_KEY = 'Overlap'

        # Conversion of GFA attributes.
        # Missing attributes types: B, J
        self.GFA_ATTRIBUTE_TYPE = {
            'i': lambda x: int(float(x)),
            'f': lambda x: float(x),
            'Z': lambda x: str(x),
            'A': lambda x: str(x),
            'H': lambda x: bytes(x)
        }
      


    #to run the tool from the file plaspipe.py
    def run(self):

        # Run the binning tool : need to see more about the parameters
        self.resultat_of_binning = self.run_binning()
        



    #the function to run the binning tool
    def run_binning(self):
        
        """ 
        args:
        command line: usually python script_binning.py
        parameters: the tool parameter
        input_bin_converted: the input file for the binning tool
        output_bin: the output file result of running the binning tool

        return: 
        resultat the binning tool
        """ 


        #extract the input format
        bin_format = self.method_configuration['input_format']  # bin_format for the format of the input of the binning tool
        
        #extract informations about the binning tool
        bin_tool_name = self.method_configuration['name']
        
        version = self.method_configuration['version']

        #extarct tool format of the output: need to find if i can do it without
        output_format = self.method_configuration['output_format']

        # Absolute path to the tool directory
        #tool_directory = 'C:/Users/user/Desktop/Bioinformatic_pipeline/plASgraph'

        # Change the working directory to the tool directory
        #os.chdir(tool_directory)

        #Formatting the input file of the pipeline 
        self.input_bin_converted = self.conversion(bin_format, bin_tool_name) #the path to input file
        print(f'here is my converted file for the binning tool {self.input_bin_converted}')
        

        # getting the tool parameter using the yaml file
        tool_parameters = self.method_configuration['parameters'] # the parameter of the binning tool

        #the file that the output results of the binning tool will be saved, have the format of the binning_tool
        output_binning = os.path.join(self.binning_dir, self.prefix + '_' + bin_tool_name + '_' + version + '.' + output_format )
        


        # Construct the full command with input_file_path and output_file
        full_command = get_command(self.method_configuration, self.input_bin_converted, output_binning)
        print(f'this is my command for the binning tool {full_command}')

        # Run the command using subprocess
        subprocess.run(full_command, shell=False)

        #need to return the output_classif
        print(f'this my binning tool result {output_binning}')
        return output_binning


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
        bin_tool_name = self.method_configuration['name']
        #prefix = self.method_configuration['prefix']
        version = self.method_configuration['version']
        


        #specify csv file path
        csv_file = csv_file = os.path.join(self.binning_dir, self.prefix + '_' + bin_tool_name + '_' + version + '.csv')
        print(f'this is my csv file enjoy {csv_file}')


        # Run the conversion script with run conversion
        run_conversion(bin_tool_name, version, file_path, csv_file)

        return csv_file

    #getter for the classification_wrapper
    def get_csv_file_from_binning(self):
        
        """ 
        returns: the csv file for the pipeline
        """

        output_pipeline_csv = self.pipeline_conversion_to_csv(self.resultat_of_binning)
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

        # If binning_dir is None, set it to the current working directory
        if self.binning_dir is None:
            self.binning_dir = os.getcwd()
        
        output_bin_file = os.path.join (self.binning_dir , self.prefix + '_converted.fasta')
        print(f"hello i wanna be converted here {output_bin_file}")

        # Create the output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_bin_file), exist_ok=True)
        

        sep = ' '
        if input_format == 'fasta':
            if self.fasta_path:
                return self.fasta_path
            elif self.gfa_path:
                print("wanna convert gfa with me?")
                #problem with the the tool name 
                
                write_GFA_to_FASTA(self.gfa_path, output_bin_file, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                print(f"let's convert this gfa file {output_bin_file}")
                return output_bin_file
                
            else:
                raise ValueError("Neither GFA nor FASTA file provided as input.")
        elif input_format == 'gfa' or input_format == 'GFA':
            if self.gfa_path:
                return self.gfa_path
            else:
                raise ValueError("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
        else:
            raise ValueError(f"You have to provide a fasta or GFA, Unsupported input format: {input_format}")
        



