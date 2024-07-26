#!/usr/bin/env python3
"""
This module contains the ClassificationWrapper class, which is responsible for managing the classification process
in a bioinformatics pipeline. It handles input formatting, running classification tools, and converting
output to the required format.

Classes:
    ClassificationWrapper: Manages the binning process, including input conversion, tool execution, and output formatting.

Dependencies:
    - os
    - logging
    - subprocess
    - shlex
    - Custom utility functions (e.g., process_exception, check_file, log_file_creation, get_command, etc.)

Usage:
    wrapper = ClassificationWrapper(classification_folder, prefix, classification_config, fasta_path, gfa_path)
    csv_output = wrapper.get_csv_file_from_classification() : standart format of the pipeline
"""
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging
import sys
import shlex
import time
import csv

from .plaspipe_utils import process_exception, process_error, create_directory, check_file, log_file_creation
from .plaspipe_utils import fasta_filter_length

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
            classification_folder (str): Path to the classification output folder.
            prefix (str): Prefix for pipeline tools.
            method_config (dict): Configuration for the classification tool.
            fasta_path (str): Path to the input FASTA file.
            gfa_path (str): Path to the input GFA file.
            gzipped_gfa (bool): Whether the GFA file is gzipped.
            gzipped_fasta (bool): Whether the FASTA file is gzipped.
        """
        self.classification_dir = classification_folder
        self.prefix = prefix
        self.method_configuration = method_config
        self.fasta_path = fasta_path
        self.gfa_path = gfa_path
        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta
        self.logger = logging.getLogger(__name__)
        self.class_tool_name = self.method_configuration['name']
        self.version = self.method_configuration['version']
        self.min_length = self.method_configuration['min_length']
        self.class_format = self.method_configuration['input_format']

    def run(self):
        """
        Run the classification tool.

        Returns:
            str: Path to the classification tool output file.
        """
        try:
            return self.run_classification()
        except Exception as e:
            process_exception(f"Error running the classification tool: {str(e)}")



    def run_classification(self):
        """
        Run the classification tool with the appropriate input and parameters.

        Returns:
            output_classification (str): Path to the classification tool output file.
        """
        try:
            # Convert input file to the required format(FASTA/GFA)
            self.input_classif_converted = self.conversion(self.class_format, self.class_tool_name)
            check_file(self.input_classif_converted)
            self.logger.info(f'Converted input file: {self.input_classif_converted}')

            output_classification = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}")

            # Exception for PlasForest tool
            if self.class_tool_name == "PlasForest" and self.version == "1.4.0":
                output_classification = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}_1.csv")
                self.logger.info(f'Output classification file for GFA input: {output_classification}')

            # Create new FASTA file containing only contigs with length > min_length
            if self.class_format == 'fasta':

                min_length = self.min_length
                if min_length == None:
                    min_length = 100

                fasta_output = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}_1.fasta")
                input_classification_file = fasta_filter_length(self.input_classif_converted, fasta_output, min_length)
                self.logger.info(f'Output classification file for FASTA input: {output_classification}')
            else:

                input_classification_file = self.input_classif_converted

            # Get the command to run the classification tool
            command = get_command(self.method_configuration, input_classification_file, output_classification)

            # Split the command string into a list if it's not already a list
            if isinstance(command, str):
                command = shlex.split(command)

            self.logger.info(f"Executing command: {' '.join(command)}")

            # Run the classification command
            result = subprocess.run(command, capture_output=False, text=True, check=True, shell=False)
            self.logger.info(f"Classification command output: {result.stdout}")

            if result.stderr:
                self.logger.warning(f"Classification command stderr: {result.stderr}")

            log_file_creation('classification_tool_result', output_classification)
            return output_classification

        except subprocess.CalledProcessError as e:
            process_error(f"Error running classification command: {e.stderr}")
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
            csv_file = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}.csv")
            run_conversion(self.class_tool_name, self.version, file_path, csv_file)
            self.logger.info(f'CSV file created: {csv_file}')
            log_file_creation('classification_pipeline_result_file', csv_file)
            return csv_file
        except Exception as e:
            process_exception(f"Error converting to CSV: {str(e)}")



    def filter_length_contig(self, file_path):
        """
        Change the contig score for contigs with length < min_length.

        Args:
            file_path (str): Path to the input CSV file.

        Returns:
            pipeline_file (str): Path to the filtered CSV file.
        """
        try:
            input_file = self.pipeline_conversion_to_csv(file_path)
            pipeline_file = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}_formated.csv")
            
            min_length = self.min_length
            if min_length == None:
                min_length = 100
                
            with open(input_file, 'r', newline='') as in_file, open(pipeline_file, 'w', newline='') as out_file:
                reader = csv.DictReader(in_file)
                fieldnames = reader.fieldnames
                writer = csv.DictWriter(out_file, fieldnames=fieldnames)
                writer.writeheader()

                for row in reader:
                    length = int(row['length'])
                    if length <= min_length:
                        row['plasmid_score'] = '0.5'
                        row['chromosome_score'] = '0.5'
                    writer.writerow(row)

            #check the new csv file
            check_file(pipeline_file)
            return pipeline_file
        except Exception as e:
            process_exception(f"Error in filter_length_contig: {str(e)}")



    def add_missing_contigs(self, file_path):
        """
        Add missing contigs (length < min_length) to the output file.

        Args:
            file_path (str): Result file of the classification tool.

        Returns:
        pipeline_file (str): Path to the updated CSV file containing all contigs, with contig scores for added contigs set to (0.5, 0.5).

        """
        try:
            input_file = self.pipeline_conversion_to_csv(file_path)
            pipeline_file = os.path.join(self.classification_dir, f"{self.prefix}_{self.class_tool_name}_{self.version}_formated.csv")
            existing_contigs = set()
            csv_contigs = []  # List to add each existing row

            # Read the CSV file
            with open(input_file, 'r', newline='') as in_file:
                reader = csv.DictReader(in_file)
                fieldnames = reader.fieldnames
                if not fieldnames:
                    raise ValueError("CSV file is empty or improperly formatted")
            
                for row in reader:
                    # Get the existing contigs from the CSV file
                    contig_name = row.get('contig_name')
                    if not contig_name:
                        raise ValueError("Missing 'contig_name' in CSV row")
                    existing_contigs.add(contig_name)
                    csv_contigs.append(row)

            # Read the FASTA file
            if not os.path.exists(self.input_classif_converted):
                raise FileNotFoundError(f"FASTA file not found: {self.input_classif_converted}")

            for record in SeqIO.parse(self.input_classif_converted, "fasta"):
                contig_name = record.id.split()[0]
                if contig_name not in existing_contigs:
                    new_contig = {
                        'contig_name': contig_name,
                        'plasmid_score': '0.5',
                        'length': str(len(record.seq)),  # Use actual sequence length
                        'chromosome_score': '0.5',
                        'label': 'ambiguous'
                        
                    }
                    csv_contigs.append(new_contig)

            # Write the updated data to pipeline file
            with open(pipeline_file, 'w', newline='') as out_file:
                writer = csv.DictWriter(out_file, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_contigs)

            check_file(pipeline_file)
            self.logger.info(f"Updated CSV file created with missing contigs: {pipeline_file}")
            return pipeline_file

        except FileNotFoundError as e:
            self.logger.error(f"File not found: {str(e)}")
            raise
        except ValueError as e:
            self.logger.error(f"Value error in CSV processing: {str(e)}")
            raise
        except IOError as e:
            self.logger.error(f"IO error when reading or writing files: {str(e)}")
            raise
        except Exception as e:
            self.logger.error(f"Unexpected error in add_missing_contigs: {str(e)}")
            raise

    def get_csv_file(self):
        """ 
        Get the CSV file for the pipeline.
        This method runs the classification, processes the results, and returns
        the path to the final CSV file, handling different input formats.

        Returns:
            str: Path to the final CSV file.
        """ 

        try:
            resultat_of_classification = self.run()
            check_file(resultat_of_classification)
            self.logger.info(f"Classification successeful. Result file: {resultat_of_classification}")

            if self.class_format == 'fasta':
                self.logger.info("Processing FASTA input format")
                pipeline_csv_file = self.add_missing_contigs(resultat_of_classification)

            elif self.class_format == 'gfa':
                self.logger.info("Processing GFA input format")
                pipeline_csv_file = self.filter_length_contig(resultat_of_classification)

            else:
                raise ValueError(f"Unsupported classification format: {self.class_format}")

            return pipeline_csv_file


        except Exception as e:
            process_exception(f"Error creating csv file: {str(e)}")


    def conversion(self, input_format, folder_tool):
        """
        Convert GFA to FASTA or return the path to the appropriate input file.

        Args:
            input_format (str): Desired input format ('fasta' or 'gfa').
            folder_tool (str): Name of the tool folder.

        Returns:
            str: Path to the input file (converted if necessary).
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
                    raise ValueError("Neither GFA nor FASTA file provided as input.")
            elif input_format.lower() == 'gfa':
                if self.gfa_path:
                    return self.gfa_path
                else:
                    raise ValueError("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
            else:
                raise ValueError(f"Unsupported input format: {input_format}")
        except Exception as e:
            process_exception(f"Error in conversion the file: {str(e)}")