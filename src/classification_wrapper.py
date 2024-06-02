import os
import subprocess
from Bio import SeqIO
from tool_conversion.classify_conversion import FastaToCsv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import sys


class ClassificationWrapper:

    """
    the classification wrapper is reponsible for the formating the input of the pipeline (fasta or GFA) and for
    converting the output of the tool format to the pipeline output format (csv)

    Args: 
    Returns:
    """ 

#initialisation
#initialisation
    def __init__(self, classification_folder, method_config, fasta_path="", gfa_path = "",gzipped_gfa=False, gzipped_fasta=False, id_fun=lambda x: x):
        #self.input_path = input_file
        #self.output_classif = output_of_classification
        self.classification_dir = classification_folder
        self.method_configuration = method_config
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


    #to run the tool 
    def run(self):


        # Run the classification tool : need to see more about the parameters
        self.resultat_of_classification = self.run_classification()
        

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
        # Extract the tool_command from the method_config dictionary
        tool_command = self.method_configuration.get("command")

        #extract the input format
        class_format = self.method_configuration['input_format']  # class_format for the format of the input of the classification tool
        
        #extract informations about the classification tool
        class_tool_name = self.method_configuration['name']
        prefix = self.method_configuration['prefix']
        version = self.method_configuration['version']

        #extarct tool format of the output: need to find if i can do it without
        output_format = self.method_configuration['output_format']


        #Formatting the input file of the pipeline 
        self.input_classif_converted = self.conversion(class_format, class_tool_name) #input_classif_converted: the path to input file
        print(f'here is my converted file {self.input_classif_converted}')
        #input_classif_converted = 'C:/Users/user/Desktop/run/input/output.fasta'

        # getting the tool parameter using the yamk file
        tool_parameters = self.method_configuration['parameters'] # the parameter of the classification tool

        #the file that the output results of the classification tool will be saved, have the format of the classification_tool
        output_classification = os.path.join(self.classification_dir, prefix + '_' + class_tool_name + '_' + version + '.' + output_format )
        

        # Split the tool_command into a list
        command_parts = tool_command.split()

        # Construct the full command with input_file_path and output_file
        full_command = command_parts + [(self.input_classif_converted), output_classification]

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
        prefix = self.method_configuration['prefix']
        version = self.method_configuration['version']

        #specify csv file path
        csv_file = csv_file = os.path.join(self.classification_dir, prefix + '_' + class_tool_name + '_' + version + '.csv')
        print(f'this is my csv file enjoy {csv_file}')
        pipeline_conv = FastaToCsv(file_path, csv_file) #path to the tool_result_file, neme of tool, version_of_tool
        pipeline_file = pipeline_conv.convert()

        return pipeline_file

    #getter for the classification_wrapper
    def get_csv_file(self):
        
        """ 
        returns: the csv file for the pipeline
        """

        output_pipeline_csv = self.pipeline_conversion_to_csv(self.resultat_of_classification)
        print("hey now I'm a csv file")
        return output_pipeline_csv


    ######################################################################
    ######################################################################
    #conversion from GFA to fasta
    def process_exception(self, msg):
        logging.exception(msg)
        print(f'EXCEPTION\t{msg}', file=sys.stderr)
        sys.exit(1)

    # Open file
    def __open_file_read(self, file_path, gzipped=False):
        """
        Open a file for reading
        """
        if gzipped:
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')

    def __open_file_write(self, file_path, gzipped=False):
        """
        Open a file for writing
        """
        if gzipped:
            return gzip.open(file_path, 'wt')
        else:
            return open(file_path, 'w')

    # Reading the fasta and gfa file
    def read_gfa(self):
        """
        Read contigs and their attributes from a GFA file
        """
        result = {}
        with self.__open_file_read(self.gfa_path, self.gzipped_gfa) as in_file:
            for gfa_line in [x for x in in_file.readlines() if x[0] == 'S']:
                line = gfa_line.rstrip()
                ctg_data = line.split('\t')
                if len(ctg_data) < 2:
                    continue  # Skip lines with fewer than 2 fields
                ctg_id = ctg_data[1]
                result[self.id_fun(ctg_id)] = ''
        return result

    def read_fasta(self):
        """ 
        args:
        fasta_path

        return:
        the contigs dict fill in with contigs name
        """ 
        print("read from fasta")
        try:
            for seq_record in SeqIO.parse(self.fasta_path, "fasta"):
                self.contigs[seq_record.id] = ""
        except Exception as e:
            self.process_exception(f'Reading {self.fasta_path}: {e}')
        return self.contigs

    # Conversion from GFA to fasta
    def __write_attributes(self, attributes_dict, keys_to_remove=[], sep=' '):
        """
        Write GFA attributes into a string

        Args:
        - attributes_dict (Dictionary): attribute key -> attribute value
        - keys_to_remove (List(str)): list of keys to not print
        - sep (str): separating string

        Returns:
        (str): list of attributes in format sep.join(key:value)
        attributes with None value are not written
        """
        return sep.join(
            [
                f'{x}:{y}'
                for x, y in attributes_dict.items()
                if x not in keys_to_remove and y is not None
            ]
        )

    def __add_attributes(self, att_data, attributes_list):
        """
        Add attributes to a dictionary

        Args:
        - att_data (List(str)): list of attribute strings in format <key>:<type>:<value>
        - attributes_list (List(str)): list of attribute keys to keep, or ['all']

        Returns:
        - (Dictionary) attribute key -> attribute value
        """
        result = {}
        for att in att_data:
            att_parts = att.split(':')
            if len(att_parts) == 3:
                key, att_type, value = att_parts
                if attributes_list == ['all'] or key in attributes_list:
                    try:
                        result[key] = self.GFA_ATTRIBUTE_TYPE.get(att_type, lambda x: x)(value)
                    except Exception as e:
                        print(f"Error processing attribute '{att}': {e}")
            else:
                print(f"Ignoring malformed attribute '{att}'")
        return result

    def read_GFA_ctgs(self, in_file_path, attributes_list, gzipped=False, ctg_fun=lambda x: x, id_fun=lambda x: x):
        """
        Read contigs and their attributes from a GFA file
        """
        result = {}
        with self.__open_file_read(in_file_path, gzipped) as in_file:
            for gfa_line in [x for x in in_file.readlines() if x[0] == 'S']:
                line = gfa_line.rstrip()
                ctg_data = line.split('\t')
                if len(ctg_data) < 2:
                    continue  # Skip lines with fewer than 2 fields
                ctg_id, ctg_seq = ctg_data[1], ctg_data[2]
                ctg_len = len(ctg_seq)
                att_data = [f'{self.GFA_SEQ_KEY}:Z:{ctg_seq}', f'{self.GFA_LEN_KEY}:i:{ctg_len}'] + ctg_data[3:]
                result[id_fun(ctg_id)] = ctg_fun(self.__add_attributes(att_data, attributes_list))
        return result

    def write_GFA_to_FASTA(self, in_GFA_file, out_FASTA_file, in_gzipped, out_gzipped, sep=' '):
        """
        Create a FASTA file from a GFA file
        """
        GFA_ctg_seqs = self.read_GFA_ctgs(in_GFA_file, attributes_list=['all'], gzipped=in_gzipped)
        ctg_records = [
            SeqRecord(
                Seq(y[self.GFA_SEQ_KEY]),
                id=x,
                name=x,
                description=f'{x}.GFA {self.__write_attributes(y, keys_to_remove=[self.GFA_SEQ_KEY])}'
            )
            for x, y in GFA_ctg_seqs.items()
        ]
        try:
            with self.__open_file_write(out_FASTA_file, gzipped=out_gzipped) as out_file:
                SeqIO.write(ctg_records, out_file, 'fasta')
        except Exception as e:
            self.process_exception(f'FASTA/GFA\tWriting {in_GFA_file} to {out_FASTA_file}: {e}')

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
        
        prefix = self.method_configuration['prefix']
        output_classif_file = os.path.join (self.classification_dir , prefix + '_converted.fasta')
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
                self.write_GFA_to_FASTA(self.gfa_path, output_classif_file, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
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
        

