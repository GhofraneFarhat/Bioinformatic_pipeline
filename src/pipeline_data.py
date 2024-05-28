from Bio import SeqIO
import gzip
import logging
import sys
import os

import shutil
from collections import defaultdict


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class PipelineData:
    def __init__(self, gfa_path, fasta_path="", output_file = "resultat.fasta", gzipped_gfa=False, gzipped_fasta=False, id_fun=lambda x: x):
        self.gfa_path = gfa_path
        self.fasta_path = fasta_path
        self.output_file = output_file
        #self.method_config = method_config
        self.gzipped_gfa = gzipped_gfa
        self.gzipped_fasta = gzipped_fasta
        self.id_fun = id_fun
        
        self.contigs = {}
        self.bins = {}
        self.others = {}  # for other types of files

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

    # Load contigs
    def load_contigs(self):
        """load contigs to the contigs dict """
        """ 
        args:
        gfa_path: the path to the file gfa input of the pipeline if it exist
        fasta_path: the path to the file fasta input of the pipeline if it exist
        """

        # GFA file
        if self.gfa_path:
            try:
                self.contigs = self.read_gfa()
                print("gfa")

            except FileNotFoundError:
                print(f"Error: GFA file '{self.gfa_path}' not found.")
            except Exception as e:
                print(f"Error reading GFA file: {e}")

        # FASTA file
        if self.fasta_path:
            try:
                self.contigs = self.read_fasta()
                print("let's use fasta")

            except FileNotFoundError:
                print(f"Error: FASTA file '{self.fasta_path}' not found.")
            except Exception as e:
                print(f"Error reading FASTA file: {e}")

    # Setter and getter of the pipelinedata
    # Update the dict of contigs after the classification method (classification wrapper)
    def set_contigs(self, contig_name, contig_class):
        self.contigs[contig_name] = {
            'contig_class': contig_class  # the score of the contig its a dict of plasmid score and chromosome score
        }

    def process_exception(self, msg):
        logging.exception(msg)
        print(f'EXCEPTION\t{msg}', file=sys.stderr)
        sys.exit(1)

    # Update the dict of bins after the binning method (binning wrapper)
    def set_bins(self, bin_id, contigs):
        self.bins[bin_id] = contigs  # contigs is a list of contig in this bin

    # Access to contigs dict
    def get_contigs(self):
        """
        return the contigs dict
        """
        return self.contigs

    # Access to the bins dict
    def get_bins(self):
        """
        return bins dict
        """
        return self.bins

    # Access to the list of contigs of a bin
    def get_bin_contig(self, bin_id):
        """
        return the list of contigs that are in bin_id
        """ 
        if bin_id in self.bins:
            return self.bins[bin_id]
        else:
            return None

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

    # The function that will provide the path to the wrapper
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
        self.output_file 
        working_dir_output = os.path.join(os.getcwd(), self.output_file)
        classify_dir_output = os.path.join(os.getcwd() + '\\' + folder_tool , self.output_file)

        sep = ' '
        if input_format == 'fasta':
            if self.fasta_path:
                return self.fasta_path
            elif self.gfa_path:
                #problem with the the tool name 
                self.write_GFA_to_FASTA(self.gfa_path, working_dir_output, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                self.write_GFA_to_FASTA(self.gfa_path, classify_dir_output, self.gzipped_gfa, self.gzipped_fasta, sep=sep)
                print("let's convert this gfa file")
                return self.output_file
                
            else:
                raise ValueError("Neither GFA nor FASTA file provided as input.")
        elif input_format == 'gfa' or input_format == 'GFA':
            if self.gfa_path:
                return self.gfa_path
            else:
                raise ValueError("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
        else:
            raise ValueError(f"You have to provide a fasta or GFA, Unsupported input format: {input_format}")
        