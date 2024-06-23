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

from .tool_command import get_command
from .tool_conversion import run_conversion
from .gfa_to_fasta import write_GFA_to_FASTA

class BinningWrapper:
    """
    The BinningWrapper is responsible for formatting the input of the pipeline (FASTA or GFA)
    and converting the output of the tool format to the pipeline output format (CSV).
    """

    def __init__(self, binning_folder, prefix, binning_config, fasta_path="", gfa_path="",  gzipped_gfa=False, gzipped_fasta=False, classification_result_file = ""):
        
        """
        Initialize the BinningWrapper
        Args:
            binning_folder (str): Path to the binning output folder (from the yaml file).
            prefix (str): Prefix
            method_config (dict): Configuration for the binning method (from the yaml file).
            fasta_path (str): Path to the input FASTA file
            gfa_path (str): Path to the input GFA file
            gzipped_gfa (bool): Whether the GFA file is gzipped
            gzipped_fasta (bool): Whether the FASTA file is gzipped
        """

        self.binning_dir = binning_folder
        self.prefix = prefix
        self.method_configuration = binning_config

        self.fasta_path = fasta_path
        self.gfa_path = gfa_path

        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta

        self.classification_result_csv = classification_result_file

    def run(self):

        """
        Run the binning tool
        Returns:
            str: Path to the binning tool output file.
        """
        try:
            return self.run_binning()
        except Exception as e:
            process_exception(f"Error running binning: {str(e)}")

    def run_binning(self):

        """
        Run the binning tool with the appropriate input and parameters
        return: 
        output_binning (str): path the the binning file result
        """

        try:
            bin_format = self.method_configuration['input_format']
            bin_tool_name = self.method_configuration['name']
            version = self.method_configuration['version']
            output_format = self.method_configuration['output_format']

            self.input_bin_converted = self.conversion(bin_format, bin_tool_name)

            #check the converted file
            check_file(self.input_bin_converted)
            print(f'Converted input file for binning: {self.input_bin_converted}')

            output_binning = os.path.join(self.binning_dir, f"{self.prefix}_{bin_tool_name}_{version}.{output_format}")
             

            full_command = get_command(self.method_configuration, self.input_bin_converted, output_binning, self.classification_result_csv)

            subprocess.run(full_command, shell=False, check=True)

            return output_binning

        except subprocess.CalledProcessError as e:
            process_error(f"Error running binning command: {str(e)}")
        except Exception as e:
            process_exception(f"Unexpected error in run_binning: {str(e)}")

    def pipeline_conversion_to_csv(self, file_path):

        """
        Convert the binning tool output to CSV format
        Return: 
        csv_file: the file of the binning result in the csv format
        """
        try:
            bin_tool_name = self.method_configuration['name']
            version = self.method_configuration['version']
            csv_file = os.path.join(self.binning_dir, f"{self.prefix}_{bin_tool_name}_{version}.csv")
            
            run_conversion(bin_tool_name, version, file_path, csv_file)
            print(f'CSV file created: {csv_file}')
            return csv_file
        except Exception as e:
            process_exception(f"Error converting to CSV: {str(e)}")

    def get_csv_file_from_binning(self):
        """
        Get the CSV file for the pipeline
        return: 
        output_pipeline_csv (str): path the csv file
        """
        try:
            resultat_of_binning = self.run()

            #check result of binning tool file
            check_file(resultat_of_binning)

            output_pipeline_csv = self.pipeline_conversion_to_csv(resultat_of_binning)
            return output_pipeline_csv

        except Exception as e:
            process_exception(f"Error getting CSV file from binning result: {str(e)}")

    def conversion(self, input_format, folder_tool):

        """
        Convert GFA to FASTA or return the path to the appropriate input file
        Return: 
        output_bin_file (str): the path to the fasta or gfa path
        """
        try:
            if self.binning_dir is None:
                self.binning_dir = os.getcwd()

            output_bin_file = os.path.join(self.binning_dir, f"{self.prefix}_converted.fasta")
            create_directory(os.path.dirname(output_bin_file), self.prefix)

            if input_format.lower() == 'fasta':
                if self.fasta_path:
                    return self.fasta_path
                elif self.gfa_path:
                    write_GFA_to_FASTA(self.gfa_path, output_bin_file, self.gzipped_gfa, self.gzipped_fasta)
                    return output_bin_file
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
            process_exception(f"Error in conversion: {str(e)}")