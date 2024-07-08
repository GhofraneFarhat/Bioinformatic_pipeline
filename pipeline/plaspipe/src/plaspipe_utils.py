import os
import logging
import sys
import gzip
import shutil

class CustomException(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)

def csv_to_tsv(csv_file_path):
    """
    Convert a CSV file to a TSV file

    Args:
        csv_file_path (str): Path to the input CSV file.

    Returns:
        str: Path to the output TSV file.
    """
    
    directory, tsv_filename = os.path.split(csv_file_path)

    # Create the output TSV filename by replacing the extension
    tsv_filename = os.path.splitext(tsv_filename)[0] + '.tsv'

    # Join the directory and TSV filename to get the output TSV file path
    tsv_file_path = os.path.join(directory, tsv_filename)

    # Convert the CSV file to TSV
    with open(csv_file_path, 'r') as csv_file, open(tsv_file_path, 'w', newline='') as tsv_file:
        next(csv_file)
        for line in csv_file:
            tsv_file.write(line.replace(',', '\t'))

    return tsv_file_path

def process_exception(msg):
    logging.exception(msg)
    print(f'EXCEPTION\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_error(msg):
    logging.error(msg)
    print(f'ERROR\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_warning(msg):
    logging.warning(msg)
    print(f'WARNING\t{msg}', file=sys.stderr)


def _check_file(in_file, log=False, msg='FILE'):
    try:
        if not os.path.isfile(in_file):
            raise CustomException('File is missing')
        elif os.path.getsize(in_file) == 0:
            process_warning(f'{msg}\t{in_file}: is empty')
    except Exception as e:
        process_exception(f'{msg}\t{in_file}: {e}')
    else:
        if log:
            logging.info(f'{msg}\t{in_file}')
            
def check_file(in_file):
    _check_file(in_file, log=False)


#gunzzipping a gfa file
def gunzip_GFA(in_file_path, out_file_path):
    """
    Gunzip a GFA file

    Args:
       in_file_path (str): path to input gzipped GFA file
       out_file_path (str): path to output GFA file

    Returns:
       Creates GFA file out_file_path
    """
    try:
        with gzip.open(in_file_path) as in_file, open(out_file_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
            print(f'this my gfa file gunzipped {out_file_path}')
            return out_file_path
    except Exception as e:
        process_exception(f'FASTA\tGunzipping {in_file_path} to {out_file_path}: {e}')


#gunzip a fasta file
def gunzip_FASTA(in_file_path, out_file_path):
    """
    Gunzip a FASTA file

    Args:
       in_file_path (str): path to input gzipped FASTA file
       out_file_path (str): path to output FASTA file

    Returns:
       Creates FASTA file out_file_path
    """
    records = []
    with gzip.open(in_file_path, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            records.append(record)
        with open(out_file_path, 'w') as out_file:
            SeqIO.write(records, out_file, 'fasta')



#verify the input fasta and gfa
def verify_input_file(gfa_path, fasta_path):
    """
    Verify that at least one input file is provided and is not empty

    Args:
        gfa_path (str): Path to the GFA file
        fasta_path (str): Path to the FASTA file

    Raises:
        ValueError: If no input file is provided or if all provided files are empty
    """
    if not gfa_path and not fasta_path:
        process_error("No input file provided. You need to provide at least one input file (GFA or FASTA)")

    files_to_check = []
    if gfa_path:
        files_to_check.append(("GFA", gfa_path))
    if fasta_path:
        files_to_check.append(("FASTA", fasta_path))

    non_empty_files = []
    for file_type, file_path in files_to_check:
        try:
            check_file(file_path)
            if os.path.getsize(file_path) > 0:
                non_empty_files.append(file_type)
        except Exception as e:
            logging.error(f"Error checking {file_type} file {file_path}: {e}")

    if not non_empty_files:
        process_error("All provided input files are empty. At least one file must be non-empty")

    logging.info(f"Valid input files: {', '.join(non_empty_files)}")




def create_directory(output_path, prefix):
    """
    Create a new path using the current working directory and a prefix if the output_path is None.
    Ensure the directory exists.

    Args:
        output_path (str): The original output path. If None, a new path will be generated.
        prefix (str): The prefix to use for the generated path.

    Returns:
        str: The original output path if not None, or a newly generated path.
    """
    if output_path is None:
        current_dir = "C:/Users/user/Desktop/Bioinformatic_pipeline/out"
        output_path = os.path.join(current_dir, f"{prefix}_output")
    
    try:
        os.makedirs(output_path, exist_ok=True)
    except OSError as e:
        logging.error(f"Error creating directory '{output_path}': {e}")
        raise

    return output_path


def check_gfa_input(class_input_format, bin_input_format, gfa_path):
    """
    Check if one of the input files is a GFA file and gfa_path is None.
    Raise an error if this condition is met.
    
    :param input_files: List of input file paths
    :param gfa_path: Path to the GFA file (can be None)
    """


    if (class_input_format.lower() == 'gfa' or bin_input_format.lower() == 'gfa') and gfa_path is None:
        error_message = "one of the tools need a gfa file, but gfa_file is None. Please provide a valid gfa_path"
        process_error(error_message)


def process_file(msg):
    process_error(f"FileNotFoundError: The file {file} was not found.")
    print(f'FILENOTFOUND\t{msg}', file=sys.stderr)
    sys.exit(1)

def exist_file(file):
    if not os.path.isfile(file):
        raise FileNotFoundError(f"YAML file '{file}' does not exist.")

def verif_file(file, format):
    if not file.endswith(format):
        raise ValueError(f"Invalid user file '{file}'. Expected a {format} file.")

# logger.py

import logging

# Setup logging
def setup_logging(log_dir):
    # Ensure the directory for the log file exists
    
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)

    logging.basicConfig(
        filename='plaspipes.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def log_file_creation(file_type, file_path):
    # Create the log entry
    log_entry = f"{file_type}: {file_path}"
    
    # Log the entry using logging.info
    logging.info(log_entry)


def get_output_directory(method_config):
    """Extract output directory and prefix from method_config"""
    outdir_pipeline = method_config['outdir_pipeline']
    prefix = method_config['prefix']
    return create_directory(outdir_pipeline, prefix), prefix

def get_gunzipped_paths(outdir_pipeline, prefix):
    """Generate paths for gunzipped files"""
    return (
        os.path.join(outdir_pipeline, f"{prefix}_gunzipped_gfa.gfa"),
        os.path.join(outdir_pipeline, f"{prefix}_gunzipped_fasta.fasta")
    )

def get_input_paths(method_config):
    """Extract input file paths from method_config"""
    return (
        method_config['input']['path_to_input_gfa'],
        method_config['input']['path_to_input_fasta']
    )

def log_input_files(gfa_path, fasta_path):
    """Log input file paths"""
    logging.info(f"GFA file: {gfa_path}")
    logging.info(f"FASTA file: {fasta_path}")
    log_file_creation('gfa_path', gfa_path)
    log_file_creation('fasta_path', fasta_path)

def process_gfa_file(gfa_path, gun_gfa_path):
    """Process GFA file, unzipping if necessary"""
    if not gfa_path:
        return None
    if not os.path.exists(gfa_path):
        raise FileNotFoundError(f"GFA file not found: {gfa_path}")
    if gfa_path.endswith('.gz'):
        try:
            return gunzip_GFA(gfa_path, gun_gfa_path)
        except IOError as e:
            raise IOError(f"Error unzipping GFA file {gfa_path}: {e}")
    elif not gfa_path.endswith('.gfa'):
        raise ValueError(f"Invalid GFA file format: {gfa_path}")
    return gfa_path

def process_fasta_file(fasta_path, gun_fasta_path):
    """Process FASTA file, unzipping if necessary"""
    if not fasta_path:
        return None
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if fasta_path.endswith('.gz'):
        try:
            return gunzip_FASTA(fasta_path, gun_fasta_path)
        except IOError as e:
            raise IOError(f"Error unzipping FASTA file {fasta_path}: {e}")
    elif not fasta_path.endswith('.fasta'):
        raise ValueError(f"Invalid FASTA file format: {fasta_path}")
    return fasta_path


""" Files and directories function """
        
def clean_files(files2clean):
    for in_file in files2clean:
        if os.path.isfile(in_file):
            os.remove(in_file)

def create_director(in_dir_list):
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

#return a default str
def process_attribute(attr):
    return "" if attr is None else attr

#return a default float
def process_attribute_float(attr):
    return 0 if attr is None else attr

def process_arguments(arg, default_arg):
    # verify the lenth
    arg = arg + [None] * (len(default_arg) - len(arg))
    
    # affect default values
    processed_arg = [default if a is None else a for a, default in zip(arg, default_arg)]
    
    return processed_arg

def absolute_path ():

    dirname = os.path.dirname(__file__)
    target_dir = "Bioinformatic_pipeline"

    # Split the path into parts
    path_parts = dirname.split(os.sep)

    # Find the index of the target directory
    try:
        target_index = path_parts.index(target_dir)
    except ValueError:
        # If target directory is not found, return the original path
        return dirname

    # Join the parts of the path up to and including the target directory
    base_path = os.sep.join(path_parts[:target_index+1])

    return base_path


