import csv
from collections import defaultdict


class gplasToCsv:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def convert(self):
        """
        convert a tab file into a CSV file

        Args:
        input_file (str): path to the input tab file.
        Returns:
        output_file (str): path to the output CSV file.
        """
        bin_dict = defaultdict(list)

        with open(self.input_file, 'r') as in_file:
        
            header = in_file.readline().strip().split()
        

            for line in in_file:
            
                fields = line.strip().split()
                if len(fields) >= 8:  # the result file contains 8 fields
                    number, contig_name, prob_chromosome, prob_plasmid, prediction, length, coverage, bin_value = fields[:8]
                    
                    gc_content = coverage

                    bin_dict[bin_value].append(f"{contig_name}:{gc_content}")
                 

    
        with open(self.output_file, 'w', newline='') as out_file:
            fieldnames = ['Bin', 'Contig']
            csv_writer = csv.DictWriter(out_file, fieldnames=fieldnames)
            csv_writer.writeheader()

            for bin_value, contigs in bin_dict.items():
                contig_str = ','.join(contigs)
                row = {'Bin': bin_value, 'Contig': contig_str}
                csv_writer.writerow(row)
            

        print(f"conversion completed {self.output_file}")
        return self.output_file

