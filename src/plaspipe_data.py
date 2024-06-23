from Bio import SeqIO
import gzip
import logging
import sys
import os
import csv
import shutil
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .plaspipe_utils import process_exception
from .plaspipe_utils import process_error
from .plaspipe_utils import process_warning
from .plaspipe_utils import check_file

from .classification_wrapper import ClassificationWrapper
from .binning_wrapper import BinningWrapper

class PipelineData:
    """
    PipelineData manages the flow of data through the pipeline. It initially contains only contig names,
    then gets updated with contig scores after classification, and finally with bin information after binning.
    """

    def __init__(self, method_config="", prefix="", classification_tool_name='classify', 
                 classification_tool_version='1.0.0', binning_tool_name='bin_tool', 
                 binning_tool_version='1.0.0', gfa_path="", fasta_path="", 
                 classification_dir="", binning_dir="", gzipped_gfa=False, 
                 gzipped_fasta=False, id_fun=lambda x: x):


        # Initialize instance variables
        self.method_configuration = method_config
        self.prefix = prefix
        self.tool_name = classification_tool_name
        self.tool_version = classification_tool_version
        self.bin_tool_name = binning_tool_name
        self.bin_tool_version = binning_tool_version
        self.gfa_path = gfa_path
        self.fasta_path = fasta_path
        self.classification_folder = classification_dir
        self.binning_folder = binning_dir
        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta
        self.id_fun = id_fun
        self.contigs = {}
        self.bins = {}
        self.others = {}  # for other types of files

    # Open file
    def __open_file_read(self, file_path, gzipped=False):
        """
        Open a file for reading
        """
        if gzipped:
            return gzip.open(file_path, 'rt')
        else:
            return open(file_path, 'r')



    def read_gfa(self):
        """
        Read contigs and their attributes from a GFA file
        return:
        result (dict): dict contains the contigs name
        """
        result = {}
        try:
            with self.__open_file_read(self.gfa_path, self.gzipped_gfa) as in_file:
                for gfa_line in [x for x in in_file.readlines() if x[0] == 'S']:
                    line = gfa_line.rstrip()
                    ctg_data = line.split('\t')
                    if len(ctg_data) < 2:
                        process_warning(f"Skipping line with fewer than 2 fields: {line}")
                        continue
                    ctg_id = ctg_data[1]
                    result[(self.tool_name, self.tool_version)] = {self.id_fun(ctg_id): ''}
            return result

        except FileNotFoundError:
            process_error(f"GFA file not found: {self.gfa_path}")
        except Exception as e:
            process_exception(f"Error reading GFA file: {e}")

    def read_fasta(self):
        """
        Read contigs from a FASTA file
        """
        try:
            for seq_record in SeqIO.parse(self.fasta_path, "fasta"):
                self.contigs[(self.tool_name, self.tool_version)] = {seq_record.id: ''}
            return self.contigs
        except FileNotFoundError:
            process_error(f"FASTA file not found: {self.fasta_path}")
        except Exception as e:
            process_exception(f"Error reading FASTA file: {e}")



    def load_contigs(self):
        """
        Load contigs to the contigs dict from either GFA or FASTA file
        """
        if self.gfa_path:
            self.contigs = self.read_gfa()
        elif self.fasta_path:
            self.contigs = self.read_fasta()
        else:
            process_error("No input file specified. Please provide either a GFA or FASTA file.")



    def update_plaspipe_data_from_classwrapper(self):
        """
        Update pipeline data using the classification tool output
        """
        try:
            classification_wrapper = ClassificationWrapper(self.classification_folder, self.prefix, 
                                                           self.method_configuration['classification'], 
                                                           self.fasta_path, self.gfa_path, 
                                                           self.gzipped_gfa, self.gzipped_fasta)
            self.output_class_pipeline = classification_wrapper.get_csv_file()
            
            #check the classification result in a csv format
            check_file(self.output_class_pipeline)

            with open(self.output_class_pipeline, 'r') as classification_result:
                reader = csv.DictReader(classification_result)
                if self.output_class_pipeline.endswith('.csv'):
                    self.update_plaspipe_data_from_class_csv(reader)
                else:
                    process_error("Unsupported output format from classification tool")
        except Exception as e:
            process_exception(f"Error in classification wrapper: {e}")



    def update_plaspipe_data_from_binwrapper(self):
        """
        Update pipeline data using the binning tool output
        """
        try:
            binning_wrapper = BinningWrapper(self.binning_folder, self.prefix, 
                                             self.method_configuration['binning'], 
                                             self.fasta_path, self.gfa_path, 
                                             self.gzipped_gfa, self.gzipped_fasta, 
                                             self.output_class_pipeline)
            output_bin_pipeline = binning_wrapper.get_csv_file_from_binning()
            
            #check the csv file from the binning wrapper 
            check_file(output_bin_pipeline)

            with open(output_bin_pipeline, 'r') as binning_result:
                reader = csv.DictReader(binning_result)
                if output_bin_pipeline.endswith('.csv'):
                    self.update_plaspipe_data_from_bin_csv(reader)
                else:
                    process_error("Unsupported output format from binning tool")
        except Exception as e:
            process_exception(f"Error in binning wrapper: {e}")



    def update_plaspipe_data_from_class_csv(self, classification_result):
        """
        Update pipeline data from classification CSV file
        """
        try:
            next(classification_result)  # Skip the header row
            for line in classification_result:
                contig_name = line['contig_name']
                plasmid_score = float(line['plasmid_score'])
                chromosome_score = float(line['chromosome_score'])
                contig_class = {
                    'plasmid_score': plasmid_score,
                    'chromosome_score': chromosome_score
                }
                self.set_contigs(self.tool_name, self.tool_version, contig_name, contig_class)
        
        except KeyError as e:
            process_error(f"Missing required column in classification CSV: {e}")
        except ValueError as e:
            process_error(f"Invalid score value in classification CSV: {e}")



    def update_plaspipe_data_from_bin_csv(self, binning_result):
        """
        Update pipeline data from binning CSV file
        """
        try:
            next(binning_result)  # Skip the header row
            for line in binning_result:
                contig_names = line['Contig']
                bin_id = line['Bin']
                contig_list = self.get_bin_contig(bin_id)
                if contig_list is None:
                    self.set_bins(bin_id, [contig_names])
                else:
                    contig_list.append(contig_names)
                    self.set_bins(bin_id, contig_list)

        except KeyError as e:
            process_error(f"Missing required column in binning CSV: {e}")


    #setter and getter
    def set_contigs(self, tool_name, tool_version, contig_name, contig_class):
        """
        Update the dict of contigs after the classification method
        """
        tool_key = (tool_name, tool_version)
        if tool_key not in self.contigs:
            self.contigs[tool_key] = {}
        self.contigs[tool_key][contig_name] = {'contig_class': contig_class}


    def set_bins(self, bin_id, contigs):
        """
        Update the dict of bins after the binning method
        """
        self.bins[bin_id] = contigs

    def get_contigs(self):
        """
        Return the contigs dict
        """
        return self.contigs

    def get_bins(self):
        """
        Return bins dict
        """
        return self.bins

    def get_bin_contig(self, bin_id):
        """
        Return the list of contigs that are in bin_id
        """
        return self.bins.get(bin_id, None)