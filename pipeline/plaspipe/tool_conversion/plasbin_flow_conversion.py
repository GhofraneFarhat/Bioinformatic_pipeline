import os

class TsvToCsv:
    def __init__(self, tsv_file, csv_file):
        self.input_file = tsv_file
        self.output_file = csv_file

    def convert(self):
        """
        Convert a TSV file to a CSV file with the same name and path, only changing the file extension.

        Returns:
            str: Path to the output CSV file.
        """

        # Open the input TSV file
        with open(self.input_file, 'r') as tsv_file:
            tsv_lines = tsv_file.readlines()[1:]  # Ignore the first line

        # Define the new column labels
        commun_labels = ['Bin', 'Flow', 'GC_bin', 'Contig']

        # Create a list of rows with the new labels
        rows = []
        rows.append(','.join(commun_labels))

        # Iterate through the lines and format each row
        #next(lines)
        for line in tsv_lines:
            fields = line.strip().split('\t')
            Bin = fields[0]
            Flow = fields[1]
            GC_bin = fields[2]
            Contig = '"' + fields[3] + '"'  # Enclose contig counts in double quotes
            row = ','.join([Bin, Flow, GC_bin, Contig])
            rows.append(row)

        # Write the rows to the output CSV file
        with open(self.output_file, 'w') as csv_file:
            csv_file.write('\n'.join(rows))

        return self.output_file

