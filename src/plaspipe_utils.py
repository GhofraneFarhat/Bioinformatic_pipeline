import os
import logging

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