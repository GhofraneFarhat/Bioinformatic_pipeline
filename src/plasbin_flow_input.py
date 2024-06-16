import os

def csv_to_tsv(csv_file_path):
    """
    Convert a CSV file to a TSV file with the same name and path, only changing the file extension.

    Args:
        csv_file_path (str): Path to the input CSV file.

    Returns:
        str: Path to the output TSV file.
    """
    # Get the directory and filename from the input CSV file path
    directory, filename = os.path.split(csv_file_path)

    # Create the output TSV filename by replacing the extension
    tsv_filename = os.path.splitext(filename)[0] + '.tsv'

    # Join the directory and TSV filename to get the output TSV file path
    tsv_file_path = os.path.join(directory, tsv_filename)

    # Convert the CSV file to TSV
    with open(csv_file_path, 'r') as csv_file, open(tsv_file_path, 'w', newline='') as tsv_file:
        for line in csv_file:
            tsv_file.write(line.replace(',', '\t'))

    return tsv_file_path

