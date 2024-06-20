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

def create_directory(directory_path):
    """
    Create a directory if it doesn't exist.

    Args:
        directory_path (str): Path to the directory.
    """
    try:
        os.makedirs(directory_path, exist_ok=True)
    except OSError as e:
        process_exception(f"Error creating directory '{directory_path}': {e}")
