import os
import logging
import sys
import csv
import gzip
import shutil
import subprocess
from .gfa_to_fasta import write_GFA_to_FASTA

# Custom exception class
class CustomException(Exception):
    def __init__(self, msg):
        # Call the base class constructor with the custom message
        super().__init__(msg)

# Conversion from CSV to TSV
def csv_to_tsv(csv_file_path):
    """
    Convert a CSV file to a TSV file

    Args:
        csv_file_path (str): Path to the input CSV file.

    Returns:
        str: Path to the output TSV file.
    """
    directory, tsv_filename = os.path.split(csv_file_path)
    tsv_filename = os.path.splitext(tsv_filename)[0] + '.tsv'
    tsv_file_path = os.path.join(directory, tsv_filename)

    try:
        with open(csv_file_path, 'r') as csv_file, open(tsv_file_path, 'w', newline='') as tsv_file:
            next(csv_file)  # Skip the header
            for line in csv_file:
                tsv_file.write(line.replace(',', '\t'))
    except Exception as e:
        process_exception(f"Error converting CSV to TSV: {e}")

    return tsv_file_path

# Exception handling functions
def process_exception(msg):
    """
    Log an exception and exit the program

    Args:
        msg (str): The exception message to log.
    """
    logging.exception(msg)
    print(f'EXCEPTION\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_error(msg):
    """
    Log an error and exit the program

    Args:
        msg (str): The error message to log.
    """
    logging.error(msg)
    print(f'ERROR\t{msg}', file=sys.stderr)
    sys.exit(1)

def process_warning(msg):
    """
    Log a warning message

    Args:
        msg (str): The warning message to log.
    """
    logging.warning(msg)
    print(f'WARNING\t{msg}', file=sys.stderr)

# File checking function
def _check_file(in_file, log=False, msg='FILE'):
    """
    Check if a file exists and is not empty

    Args:
        in_file (str): Path to the file.
        log (bool): Whether to log the file check.
        msg (str): Message to log.
    """
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
    """
    Check if a file exists and is not empty, without logging

    Args:
        in_file (str): Path to the file.
    """
    _check_file(in_file, log=False)

def log_file(in_file):
    """
    Check if a file exists and is not empty, with logging

    Args:
        in_file (str): Path to the file.
    """
    _check_file(in_file, log=True)

# Gunzip a GFA file
def gunzip_GFA(in_file_path, out_file_path):
    """
    Gunzip a GFA file

    Args:
        in_file_path (str): Path to input gzipped GFA file
        out_file_path (str): Path to output GFA file

    Returns:
        str: Path to the output GFA file
    """
    try:
        with gzip.open(in_file_path, 'rb') as in_file, open(out_file_path, 'wb') as out_file:
            shutil.copyfileobj(in_file, out_file)
        print(f'GFA file gunzipped to {out_file_path}')
        return out_file_path
    except Exception as e:
        process_exception(f'Error gunzipping GFA file {in_file_path} to {out_file_path}: {e}')
# Gunzip a FASTA file
def gunzip_FASTA(in_file_path, out_file_path):
    """
    Gunzip a FASTA file

    Args:
        in_file_path (str): Path to input gzipped FASTA file
        out_file_path (str): Path to output FASTA file

    Returns:
        str: Path to the output FASTA file
    """
    try:
        records = []
        with gzip.open(in_file_path, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                records.append(record)
        with open(out_file_path, 'w') as out_file:
            SeqIO.write(records, out_file, 'fasta')
        return out_file_path
    except Exception as e:
        process_exception(f'Error gunzipping FASTA file {in_file_path} to {out_file_path}: {e}')
# Verify input files
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
# Create a directory
def create_directory(output_path, prefix):
    """
    Create a new path using the current working directory and a prefix if the output_path is None. Ensure the directory exists.

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
# Check GFA input
def check_gfa_input(class_input_format, bin_input_format, gfa_path):
    """
    Check if one of the input files is a GFA file and gfa_path is None. Raise an error if this condition is met.

    Args:
        class_input_format (str): Class input format
        bin_input_format (str): Binary input format
        gfa_path (str): Path to the GFA file (can be None)
    """
    if (class_input_format.lower() == 'gfa' or bin_input_format.lower() == 'gfa') and gfa_path is None:
        error_message = "One of the tools needs a GFA file, but gfa_file is None. Please provide a valid gfa_path"
        process_error(error_message)
# File existence check
def process_file(msg):
    """
    Process a file error message

    Args:
        msg (str): The error message to process.
    """
    process_error(f"FileNotFoundError: The file {file} was not found.")
    print(f'FILENOTFOUND\t{msg}', file=sys.stderr)
    sys.exit(1)

def exist_file(file):
    """
    Check if a file exists

    Args:
        file (str): Path to the file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    if not os.path.isfile(file):
        raise FileNotFoundError(f"YAML file '{file}' does not exist.")

def verif_file(file, format):
    """
    Verify the file format

    Args:
        file (str): Path to the file.
        format (str): Expected file format.

    Raises:
        ValueError: If the file format is incorrect.
    """
    if not file.endswith(format):
        raise ValueError(f"Invalid user file '{file}'. Expected a {format} file.")       

# Setup logging
def setup_logging(log_dir):
    """
    Setup logging

    Args:
        log_dir (str): Directory for log files
    """
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logging.basicConfig(
        filename='plaspipes.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def log_file_creation(file_type, file_path):
    """
    Log file creation

    Args:
        file_type (str): Type of the file
        file_path (str): Path to the file
    """
    log_entry = f"{file_type}: {file_path}"
    logging.info(log_entry)

# Extract output directory and prefix from method_config
def get_output_directory(method_config):
    """
    Extract output directory and prefix from method_config

    Args:
        method_config (dict): Configuration dictionary

    Returns:
        tuple: Output directory and prefix
    """
    outdir_pipeline = method_config['outdir_pipeline']
    prefix = method_config['prefix']
    return create_directory(outdir_pipeline, prefix), prefix

# Generate paths for gunzipped files
def get_gunzipped_paths(outdir_pipeline, prefix):
    """
    Generate paths for gunzipped files

    Args:
        outdir_pipeline (str): Output directory for the pipeline
        prefix (str): Prefix for the files

    Returns:
        tuple: Paths for gunzipped GFA and FASTA files
    """
    return (
        os.path.join(outdir_pipeline, f"{prefix}_gunzipped_gfa.gfa"),
        os.path.join(outdir_pipeline, f"{prefix}_gunzipped_fasta.fasta")
    )

# Extract input file paths from method_config
def get_input_paths(method_config):
    """
    Extract input file paths from method_config

    Args:
        method_config (dict): Configuration dictionary

    Returns:
        tuple: Paths to input GFA and FASTA files
    """
    return (
        method_config['input']['path_to_input_gfa'],
        method_config['input']['path_to_input_fasta']
    )

# Log input file paths
def log_input_files(gfa_path, fasta_path):
    """
    Log input file paths

    Args:
        gfa_path (str): Path to the GFA file
        fasta_path (str): Path to the FASTA file
    """
    logging.info(f"GFA file: {gfa_path}")
    logging.info(f"FASTA file: {fasta_path}")
    log_file_creation('gfa_path', gfa_path)
    log_file_creation('fasta_path', fasta_path)

# Process and verify GFA file as pipeline input, unzipping if necessary
def process_gfa_file(gfa_path, gun_gfa_path):
    """
    Process and verify GFA file as pipeline input, unzipping if necessary

    Args:
        gfa_path (str): Path to the GFA file
        gun_gfa_path (str): Path to the gunzipped GFA file

    Returns:
        str: Path to the processed GFA file
    """
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

# Process and verify the FASTA file, unzipping if necessary
def process_fasta_file(fasta_path, gun_fasta_path):
    """
    Process and verify the FASTA file, unzipping if necessary

    Args:
        fasta_path (str): Path to the FASTA file
        gun_fasta_path (str): Path to the gunzipped FASTA file

    Returns:
        str: Path to the processed FASTA file
    """
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

# Clean files
def clean_files(files2clean):
    """
    Remove specified files

    Args:
        files2clean (list): List of file paths to remove
    """
    for in_file in files2clean:
        if os.path.isfile(in_file):
            os.remove(in_file)

# Create directories
def create_director(in_dir_list):
    """
    Create directories if they do not exist

    Args:
        in_dir_list (list): List of directory paths to create
    """
    for in_dir in in_dir_list:
        if not os.path.exists(in_dir):
            os.makedirs(in_dir)

# Process attribute
def process_attribute(attr):
    """
    Return a default string if the attribute is None

    Args:
        attr: Attribute value

    Returns:
        str: Processed attribute
    """
    return "" if attr is None else attr

# Process attribute as float
def process_attribute_float(attr):
    """
    Return a default float if the attribute is None

    Args:
        attr: Attribute value

    Returns:
        float: Processed attribute
    """
    return 0 if attr is None else attr

# Process arguments
def process_arguments(arg, default_arg):
    """
    Assign default parameters

    Args:
        arg (list): List of arguments
        default_arg (list): List of default arguments

    Returns:
        list: Processed arguments
    """
    arg = arg + [None] * (len(default_arg) - len(arg))
    processed_arg = [default if a is None else a for a, default in zip(arg, default_arg)]
    return processed_arg

# Get absolute path of the project
def absolute_path():
    """
    Get the absolute path of the project

    Returns:
        str: Absolute path of the project
    """
    dirname = os.path.dirname(__file__)
    target_dir = "Bioinformatic_pipeline"
    path_parts = dirname.split(os.sep)
    try:
        target_index = path_parts.index(target_dir)
    except ValueError:
        return dirname
    base_path = os.sep.join(path_parts[:target_index+1])
    return base_path

# Run external command
def _run_cmd(cmd, output, num_attempts, exit_on_error):
    """
    Run external command, trying at most num_attempts times

    Args:
        cmd (list): Command to run
        output (str): Output file path
        num_attempts (int): Number of attempts
        exit_on_error (bool): Exit on error flag

    Returns:
        int: Process return code
    """
    cmd_str = ' '.join(cmd)
    logging.info(f'COMMAND\t{cmd_str}')
    attempt = 1
    process_returncode = -1
    while attempt <= num_attempts:
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            msg = f'COMMAND\t{cmd_str} attempt #{attempt} {e}'
            if attempt < num_attempts:
                process_warning(f'{msg}: retrying')
            elif exit_on_error:
                process_exception(f'{msg}: aborting')
            else:
                process_warning(f'{msg}: failed but not aborting')
        else:
            if output is None:
                logging.info(f'STDOUT:\n{process.stdout}')
            else:
                with open(output, 'w') as out_file:
                    out_file.write(process.stdout)
            if len(process.stderr) > 0:
                logging.warning(f'STDERR:\n{process.stderr}')
            process_returncode = process.returncode
        attempt += 1
    return process_returncode

def run_cmd(cmd, num_attempts=5, exit_on_error=True):
    """
    Run external command, trying at most num_attempts=5 times

    Args:
        cmd (list): Command to run
        num_attempts (int): Number of attempts
        exit_on_error (bool): Exit on error flag

    Returns:
        int: Process return code
    """
    return _run_cmd(cmd, None, num_attempts, exit_on_error)

# Convert GFA to FASTA or return the path to the appropriate input file
def conversion_gfa_fasta(gfa_file, fasta_file, gzipped_gfa = False, gzipped_fasta = False):
    """
    Convert GFA to FASTA or return the path to the appropriate input file.

    Args:
        input_format (str): Desired input format ('fasta' or 'gfa').
        folder_tool (str): Name of the tool folder.

    Returns:
        str: Path to the input file (converted if necessary).
    """
    try:
        sep = ' '
        write_GFA_to_FASTA(gfa_file, fasta_file, gzipped_gfa, gzipped_fasta, sep=sep)
        log_file_creation('gfa converted file', fasta_file)
        return fasta_file

    except Exception as e:
        process_exception(f"Error in conversion the file: {str(e)}")
#copy file to folder for rfplasmid
def copy_file(source_file, destination_folder):
    """ 
    copy file into a specific destination
    """ 

    try:
        shutil.copy(source_file, destination_folder)
        print(f"File '{source_file}' copied to '{destination_folder}' successfully.")
    except IOError as error:
        print(f"Error to copy file: {error}")
#filter the short contigs for fasta file
def fasta_filter_length(input_file, output_file, min_leng):


    """ 
    create a new fasta file contain only contigs with length > min_leng 
    Args:
    input_file: the fasta file (user input/ gfa converted)
    min_leng: int provodided by the user
    Return:
    output_file: fasta file contains contigs > min_leng
    """
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        current_header = ''
        current_sequence = ''
        
        for line in in_file:
            line = line.strip()
            if line.startswith('>'):
                # Process the previous contig if it exists
                if current_header and len(current_sequence) > min_leng:
                    out_file.write(f"{current_header}\n{current_sequence}\n")
                
                # Start a new contig
                current_header = line
                current_sequence = ''
            else:
                current_sequence += line
        
        # Process the last contig
        if current_header and len(current_sequence) > min_leng:
            out_file.write(f"{current_header}\n{current_sequence}\n")

    return output_file
#conversion from csv to tab
def csv_to_tab(input_csv):
    """
    Convert a CSV file to a TAB file

    Args:
        input_csv (str): Path to the input CSV file.

    Returns:
        output_tab (str): Path to the output TSV file.
    """

    directory, tab_filename = os.path.split(input_csv)

    # Create the output TSV filename by replacing the extension
    tab_filename = os.path.splitext(tab_filename)[0] + '_1' + '.tab'

    # Join the directory and TSV filename to get the output TSV file path
    output_tab = os.path.join(directory, tab_filename)

    with open(input_csv, 'r') as in_file, open(output_tab, 'w') as out_file:
        # Write the header for the output tab file
        out_file.write("Prob_Chromosome\tProb_Plasmid\tPrediction\tContig_name\tContig_length\n")
        
        # Read the input CSV file
        csv_reader = csv.DictReader(in_file)
        
        for row in csv_reader:
            contig_name = f"{row['contig_name']}"
            chromosome_score = float(row['chromosome_score'])
            plasmid_score = float(row['plasmid_score'])
            length = int(row['length'])
            
            # Determine the prediction
            prediction = "Chromosome" if chromosome_score > plasmid_score else "Plasmid"
            
            # Write the transformed data to the output file
            out_file.write(f"{chromosome_score:.6f}\t{plasmid_score:.6f}\t{prediction}\t{contig_name}\t{length}\n")

    return output_tab
#add eliminated short contigs
def update_contig_names(input_tab, gfa_file):
    """
    Updates the contig names in the input tab file based on the information provided in the GFA file.

    Parameters:
    - input_tab: str, path to the input tab file.
    - gfa_file: str, path to the GFA file containing new contig information.
    - output_tab: str, path to the output tab file with updated contig names.

    Returns:
    - tab file
    """

    directory, tab_filename = os.path.split(input_tab)

    # Create the output TSV filename by replacing the extension
    tab_filename = 'gfa_names_gplas2.tab'

    # Join the directory and TSV filename to get the output TSV file path
    output_tab = os.path.join(directory, tab_filename)

    # Read the GFA file and create a mapping of old contig names to new contig names
    contig_mapping = {}
    
    with open(gfa_file, 'r') as gfa:
        for line in gfa:
            if line.startswith('S'):
                parts = line.split()
                old_contig_name = parts[1]  # The old contig name (number)
                length = parts[3].split(':')[2]  # Extract length from LN:i:114
                depth = parts[4].split(':')[2]  # Extract depth from dp:f:2.789286085510359
                new_contig_name = f"S{old_contig_name}_LN:i:{length}_dp:f:{depth}"
                contig_mapping[old_contig_name] = new_contig_name

    # Create the output tab file with updated contig names
    with open(input_tab, 'r') as infile, open(output_tab, 'w') as outfile:
        # Write the header for the output file
        outfile.write("Prob_Chromosome\tProb_Plasmid\tPrediction\tContig_name\tContig_length\n")
        
        # Read the input tab file
        tab_reader = csv.DictReader(infile, delimiter='\t')
        
        for row in tab_reader:
            old_contig_name = row['Contig_name']
            # Update the contig name using the mapping
            new_contig_name = contig_mapping.get(old_contig_name, old_contig_name)  # Fallback to old name if not found
            
            # Write the updated row to the output file
            outfile.write(f"{row['Prob_Chromosome']}\t{row['Prob_Plasmid']}\t{row['Prediction']}\t{new_contig_name}\t{row['Contig_length']}\n")

    return output_tab
#gzip files
def gzip_file(file_path, binning_outdir = ""):
    """
    gzip the given file and creates a new file with '.gz' extension
    
    Args:
        file_path (str): Path to the file to be gzipped
    """
    # Get the file name and folder path
    directory, file_name = os.path.split(file_path)
    
    # Create the gzipped file path
    gzipped_file_path = os.path.join(binning_outdir, file_name + '.gz')
    
    # Open the input file in binary mode
    with open(file_path, 'rb') as input_file:
        # Create the gzipped file and write the compressed data
        with gzip.open(gzipped_file_path, 'wb') as gzipped_file:
            shutil.copyfileobj(input_file, gzipped_file)

    return gzipped_file_path
    
    print(f"File '{file_path}' has been gzipped to '{gzipped_file_path}'")