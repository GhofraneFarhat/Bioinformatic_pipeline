import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import sys

from .plaspipe_utils import process_exception
from .plaspipe_utils import process_error
from .plaspipe_utils import create_directory
from .plaspipe_utils import check_file
from .plaspipe_utils import log_file_creation

from .tool_command import get_command
from .tool_conversion import run_conversion
from .gfa_to_fasta import write_GFA_to_FASTA

class ClassificationWrapper:
    """
    The ClassificationWrapper is responsible for formatting the input of the pipeline (FASTA or GFA)
    and converting the output of the tool format to the pipeline output format (CSV).
    """

    def __init__(self, classification_folder, prefix, method_config, fasta_path="", gfa_path="", gzipped_gfa=False, gzipped_fasta=False):
        """
        Initialize the ClassificationWrapper.

        Args:
            classification_folder (str): Path to the classification output folder (from the yaml file).
            prefix (str): Prefix
            method_config (dict): Configuration for the classification method (from the yaml file).
            fasta_path (str): Path to the input FASTA file
            gfa_path (str): Path to the input GFA file
            gzipped_gfa (bool): Whether the GFA file is gzipped
            gzipped_fasta (bool): Whether the FASTA file is gzipped
        """
        self.classification_dir = classification_folder
        self.prefix = prefix
        self.method_configuration = method_config
        self.fasta_path = fasta_path
        self.gfa_path = gfa_path
        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta

    def run(self):
        """
        Run the classification tool.

        Returns:
            str: Path to the classification tool output file.
        """
        try:
            result_of_wrapper = self.run_classification()
            return result_of_wrapper
        except Exception as e:
            process_exception(f"Error running the classification tool: {str(e)}")

    def run_classification(self):
        """
        Run the classification tool with the appropriate input and parameters.

        Returns:
            output_classification (str): Path to the classification tool output file.
        """
        try:
            class_format = self.method_configuration['input_format']
            class_tool_name = self.method_configuration['name']
            version = self.method_configuration['version']
            

            #need to fix it
            tool_directory = 'C:/Users/user/Desktop/Bioinformatic_pipeline/plASgraph'
            os.chdir(tool_directory)

            self.input_classif_converted = self.conversion(class_format, class_tool_name)

            #check the converted file 
            check_file(self.input_classif_converted)
            logging.info(f'Converted input file: {self.input_classif_converted}')
            print(f'Converted input file: {self.input_classif_converted}')

            
            output_classification = os.path.join(self.classification_dir, f"{self.prefix}_{class_tool_name}_{version}")

            full_command = get_command(self.method_configuration, self.input_classif_converted, output_classification)

            subprocess.run(full_command, shell=False, check=True)

            log_file_creation('classifcation_tool_result', output_classification)
            return output_classification

        except subprocess.CalledProcessError as e:
            process_error(f"Error running classification command: {str(e)}")
        except Exception as e:
            process_exception(f"Unexpected error in run_classification: {str(e)}")

    def pipeline_conversion_to_csv(self, file_path):
        """
        Convert the classification tool output to CSV format.

        Args:
            file_path (str): Path to the classification tool output file.

        Returns:
            csv_file (str): Path to the converted CSV file.
        """
        try:
            class_tool_name = self.method_configuration['name']
            version = self.method_configuration['version']
            csv_file = os.path.join(self.classification_dir, f"{self.prefix}_{class_tool_name}_{version}.csv")
            
            run_conversion(class_tool_name, version, file_path, csv_file)
            print(f'CSV file created: {csv_file}')
            log_file_creation('classifcation_pipeline_result', csv_file)
            return csv_file
        except Exception as e:
            process_exception(f"Error converting to CSV: {str(e)}")

    def get_csv_file(self):
        """
        Get the CSV file for the pipeline.

        Returns:
            output_pipeline_csv (str): Path to the CSV file.
        """
        try:
            resultat_of_classification = self.run()

            #chech classifcation result file
            check_file(resultat_of_classification)

            output_pipeline_csv = self.pipeline_conversion_to_csv(resultat_of_classification)
            return output_pipeline_csv

        except Exception as e:
            process_exception(f"Error getting CSV file: {str(e)}")

    def conversion(self, input_format, folder_tool):
        """
        Convert GFA to FASTA or return the path to the appropriate input file.

        Args:
            input_format (str): Desired input format ('fasta' or 'gfa').
            folder_tool (str): Name of the tool folder.

        Returns:
            output_classif_file (str): Path to the input file (converted if necessary).
        """
        try:
            if self.classification_dir is None:
                self.classification_dir = os.getcwd()

            output_classif_file = os.path.join(self.classification_dir, f"{self.prefix}_converted.fasta")
            create_directory(os.path.dirname(output_classif_file), self.prefix)

            sep = ' '

            if input_format.lower() == 'fasta':
                if self.fasta_path:
                    return self.fasta_path
                elif self.gfa_path:
                    write_GFA_to_FASTA(self.gfa_path, output_classif_file, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                    log_file_creation('classification_tool_input', output_classif_file)
                    return output_classif_file
                else:
                    process_error("Neither GFA nor FASTA file provided as input.")
            elif input_format.lower() == 'gfa':
                if self.gfa_path:
                    return self.gfa_path
                else:
                    process_error("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
            else:
                process_error(f"Unsupported input format: {input_format}")
        except Exception as e:
            process_exception(f"Error in conversion the file: {str(e)}")