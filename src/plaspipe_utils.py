import os
import logging
import sys

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


