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

#from classification_wrapper import ClassificationWrapper

class PipelineData:

    """ 
    plaspipe_data is the flow for the pipeline, firstly before the wrapper, the data contains only the contigs name
    then after the classifcation wrapper it gets updated with the contigs score, after the binning wrapper 
    it get updated with the bins
    """ 
    def __init__(self, classification_wrapper,binning_wrapper, classification_tool_name, classification_tool_version, gfa_path= "", fasta_path="", gzipped_gfa=False, gzipped_fasta=False, id_fun=lambda x: x):
     
        self.classification_wrapper = classification_wrapper
        self.binning_wrapper = binning_wrapper

        self.tool_name = classification_tool_name
        self.tool_version = classification_tool_version

        self.gfa_path = gfa_path
        self.fasta_path = fasta_path
        
        #self.method_config = method_config
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
                result[(self.tool_name, self.tool_version)] = {self.id_fun(ctg_id): ''}
                
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
                self.contigs[(self.tool_name, self.tool_version)] = {seq_record.id: ''}

        except Exception as e:
            self.process_exception(f'Reading {self.fasta_path}: {e}')
        return self.contigs

    # Load contigs
    def load_contigs(self):
        
        """
        load contigs to the contigs dict
        Args:
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




    #the pipeline for updating, will use only a standart format: csv format that contains contig_name,chromosome_score, plasmid_score
    #let's say that this file called: pipeline_file_resultat.csv
    #the function to update the pipeline instance
    # update the pipeline data using the classification tool, classification_wrapper provide the plaspipe_data with a csv file 
    def update_plaspipe_data_from_classwrapper(self): #pipeline_file_resultat
        print("let's update")
        
        #convert the tool output format to the pipeline output format
        output_class_pipeline = self.classification_wrapper.get_csv_file()
        print("let's use a csv file")
        print(output_class_pipeline)


        with open(output_class_pipeline, 'r') as classification_result:
            reader = csv.DictReader(classification_result)
            print(f"I'm using {output_class_pipeline}")


            if output_class_pipeline.endswith ('.csv'): #we will need just the csv format 
                self.update_plaspipe_data_from_class_csv(reader)


            else:
                print("Unsupported output format")

    def update_plaspipe_data_from_binwrapper(self): #pipeline_file_resultat
        print("let's update")
        
        #convert the tool output format to the pipeline output format
        output_bin_pipeline = self.binning_wrapper.get_csv_file()
        print("let's use a csv file")
        print(output_bin_pipeline)


        with open(output_bin_pipeline, 'r') as binning_result:
            reader = csv.DictReader(binning_result)
            print(f"I'm using {output_bin_pipeline}")


            if output_bin_pipeline.endswith ('.csv'): #we will need just the csv format 
                self.update_plaspipe_data_from_bin_csv(reader)


            else:
                print("Unsupported output format")



    # update from a csv file
    #need to change the update from csv to use the csv biblio
    def update_plaspipe_data_from_class_csv(self, classification_result):

               
        next(classification_result)# Skip the header row if there is one

        for line in classification_result:
            
            contig_name = line['contig_name']
            plasmid_score = line['plasmid_score']
            chromosome_score = line['chromosome_score']

            contig_class = {
                'plasmid_score': float(plasmid_score),
                'chromosome_score': float(chromosome_score)
            }#dict of contigs from the classification tool

            self.set_contigs(self.tool_name, self.tool_version, contig_name, contig_class)

    def update_plaspipe_data_from_class_bin_csv(self, binning_result):

        # Skip the first line
        next(binning_result)
        
        for line in binning_result:

            contig_name = line['contig_name']
            bin_id = line['bin_id']


            contig_list = self.pipeline_data.get_bin_contig(bin_id)
            if contig_list is None:
                self.pipeline_data.set_bins(bin_id, [contig_name])
            else:
                contig_list.append(contig_name)
                self.pipeline_data.set_bins(bin_id, contig_list)

    #we need to handle if in the conversation of the file we lose a contig cause of an error 


    


    # Setter and getter of the pipelinedata
    # Update the dict of contigs after the classification method (classification wrapper)
    def set_contigs(self, tool_name, tool_version, contig_name, contig_class):

        tool_key = (tool_name, tool_version)

        if tool_key not in self.contigs:
            self.contigs[tool_key] = {}

        self.contigs[tool_key][contig_name] = {
            'contig_class': contig_class  # the score of the contig its a dict of plasmid score and chromosome score
        }



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

    #the exceptions functions
    def process_exception(self, msg):
        logging.exception(msg)
        print(f'EXCEPTION\t{msg}', file=sys.stderr)
        sys.exit(1)



        
""" 
# Example usage
fasta_file = 'output.fasta'
in_gfa_file = 'binniginput.GFA'
method_config = {"name": "classify", "parameters":"", "command":"python classify.py", "input_format": "fasta"}

# Initialize the ClassificationWrapper object
# Make sure to pass the arguments in the correct order
classy = ClassificationWrapper("res.fasta", 
                                method_config, 
                                "C:/Users/user/Desktop/run/classify", 
                                "conv.fasta", 
                                fasta_path="", 
                                gfa_path="binniginput.GFA")

classy.run()
pip = PipelineData(classy, "binniginput.GFA", "")
pip.load_contigs()
contigs = pip.get_contigs()
pip.update_plaspipe_data_from_classwrapper()

contigss = pip.get_contigs()
print(contigs)
print(contigss)

"""