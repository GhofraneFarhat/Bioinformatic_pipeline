import os

class TsvToCsv:
    def __init__(self, tsv_file, csv_file):
        self.input_file = tsv_file
        self.output_file = csv_file

    def convert(self):
        """
        convert a TSV file to a CSV file
        Args:
            tsv_file (str): path to tsv file 
        Returns:
            csv_file (str): Path to the output CSV file.
        """

        with open(self.input_file, 'r') as tsv_file:
            tsv_lines = tsv_file.readlines()[1:]  #pass header

        
        csv_labels = ['Bin', 'Flow', 'GC_bin', 'Contig']
        rows = []
        rows.append(','.join(csv_labels))

        for line in tsv_lines:
            fields = line.strip().split('\t')
            Bin = fields[0]
            Flow = fields[1]
            GC_bin = fields[2]
            Contig = '"' + fields[3] + '"'  # 
            row = ','.join([Bin, Flow, GC_bin, Contig])
            rows.append(row)

        
        with open(self.output_file, 'w') as csv_file:
            csv_file.write('\n'.join(rows))


        print(f'conversion of plasbin-flow is done {self.output_file}')
        return self.output_file

