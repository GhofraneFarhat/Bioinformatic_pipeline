from Bio import SeqIO
import gzip
import logging
import subprocess
import os


class PipelineData:
    def __init__(self, gfa_path = "", fasta_path = "", gfa_gzipped=False, fasta_gzipped=False):
        self.gfa_path = gfa_path
        self.fasta_path = fasta_path
        self.gfa_gzipped = gfa_gzipped
        self.fasta_gzipped = fasta_gzipped
        self.contigs = {}
        self.bins = {}
        self.others={} #des autre types de fichiers

#functions to check the file 
    """
    #verify if the input file (fasta or GFA file exists before passing it to the binning or the classification method)
    def input_checker(self, input_file):
        if not os.path.isfile(input_file):
            raise argparse.ArgumentTypeError(f'Invalid Argument! {input_file} does not exist.')
            print("1")
        return input_file
        


    #verify the GFA file or fasta file
    def process_exception(msg):
        print("2")
        logging.exception(msg)
        print(f'EXCEPTION\t{msg}', file=sys.stderr)
        sys.exit(1)


    def process_warning(msg):
        print("3")
        logging.warning(msg)
        print(f'WARNING\t{msg}', file=sys.stderr)

         #Files and directories function 


    #verifier l'existance du fichier
    def _check_file(in_file, log=False, msg='FILE'):
        try:
            if not os.path.isfile(in_file):
                raise CustomException('File is missing')
                print("4")
            elif os.path.getsize(in_file) == 0:
                process_warning(f'{msg}\t{in_file}: is empty')
        except Exception as e:
            process_exception(f'{msg}\t{in_file}: {e}')
        else:
            if log:
                logging.info(f'{msg}\t{in_file}')
            
    def check_file(in_file):
        _check_file(in_file, log=False)    


    # for the fasta and GFA, we can add a method to verify if the conversion of the files is well done or not
    #i need to add the gzipped file
    def verify (self, fasta_path, gfa_path):
        fasta_contigs = []
        GFA_contigs = []
        fasta_contigs = self.read_fasta(fasta_path)
        GFA_contigs = self.read_gfa(gfa_path)

        if GFA_contigs == fasta_contigs:
            print("okey")
        else:
            print("error")

    #maybe i can modifiy the print of errors with a logger 



    """

    def load_contigs(self):
        # GFA file

        if self.gfa_path:
            try:
                if self.gfa_path.endswith('.gz'):
                    with gzip.open(self.gfa_path, 'rt') as gfa_file:
                        self.read_gfa(gfa_file)
                else:
                    with open(self.gfa_path, 'r') as gfa_file:
                        self.read_gfa(gfa_file)
            except FileNotFoundError:
                print(f"Error: GFA file '{self.gfa_path}' not found.")
            except Exception as e:
                print(f"Error reading GFA file: {e}")



        # FASTA file
        if self.fasta_path:
            try:
                if self.fasta_path.endswith('.gz'):
                    with gzip.open(self.fasta_path, 'rt') as fasta_file:
                        self.read_fasta(fasta_file)
                else:
                    with open(self.fasta_path, 'r') as fasta_file:
                        self.read_fasta(fasta_file)
            except FileNotFoundError:
                print(f"Error: FASTA file '{self.fasta_path}' not found.")
            except Exception as e:
                print(f"Error reading FASTA file: {e}")

        print("contigs loaded")

#still having problems with the read GFA 
    def read_gfa(self, gfa_path):
        pass

#read the fasta file 
    def read_fasta(self, fasta_file):
        try:
            for seq_record in SeqIO.parse(fasta_file, "fasta"):
                self.contigs[seq_record.id] = ""
        except ValueError as e:
            print(f"Error parsing FASTA file: {e}")




    #update the dict of contigs after the classification method (classification wrapper)
    def set_contigs(self, contig_name, contig_class):
        self.contigs[contig_name] = {

            'contig_class': contig_class #the score of the contig its a dict of plasmid score and chromosome score

        }


    #update the dict of bins after the binning method (binning wrapper)
    def set_bins(self, bin_id, contigs):
        self.bins[bin_id] = contigs #contigs is a list of contig in this bin


    #access to contigs dict
    def get_contigs(self):
        return self.contigs

    #access to the bins dict
    def get_bins(self):
        return self.bins

    #access to the list of contigs of a bin
    def get_bin_contig(self, bin_id):
        
        if bin_id in self.bins:
            return self.bins[bin_id]
        else:
            return None


# pipelinedata.py
    def convert_gfa_to_fasta(self, gfa_path):
        fasta_out = ""
        num_seqs = 0
        with open(gfa_path, "r") as gfa_file:

            for line in gfa_file:
                if line.startswith("S\t"):
                    split = line.strip().split("\t")
                    seq = split[2]
                    fasta_out += ">{}\n".format(split[1])
                    fasta_out += split[2] + "\n"
                    num_seqs += 1

        fasta_file = os.path.splitext(gfa_path)[0] + '.fasta'

        with open(fasta_file, "w") as fasta_file:
            fasta_file.write(fasta_out)

        return (fasta_file)

#the function that will provide the path to the wrapper

    def conversion(self, input_format):
        
        if input_format == 'fasta':
            if self.fasta_path:
                return self.fasta_path
            elif self.gfa_path:
                fasta_path = self.convert_gfa_to_fasta(self.gfa_path)
                return fasta_path 
            else:
                raise ValueError("Neither GFA nor FASTA file provided as input.")
        elif input_format == 'gfa'or input_format == 'GFA':
            if self.gfa_path:
                return self.gfa_path
            else:
                raise ValueError("GFA file not provided as input, and conversion from FASTA to GFA is not possible.")
        else:
            raise ValueError(f"You have to provide a fasta or GFA, Unsupported input format: {input_format}")

