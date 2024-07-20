import csv
import shutil
import argparse

class Plasgraph2ToCsv:
    def __init__(self, csv_input, csv_output):
        
        self.input_file = csv_input
        self.output_file = csv_output

    def convert(self):
        """
        convert a csv plasgraph2 file to  CSV file, for plasgraph2 it's just row changement the tool output is already a CSV file

        Args:
        csv_file (str): path to the input csv file.
        Returns:
        output_csv (str): path to the output CSV file.
        """
        with open(self.input_file, 'r') as file:
            reader = csv.reader(file)
            data = list(reader)

        
        new_data = []
        for row in data[1:]:  
            new_row = [row[1], row[3],row[2], row[4], row[5]]
            new_data.append(new_row)

        
        with open(self.output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['contig_name','plasmid_score','length','chromosome_score','label'])  
            writer.writerows(new_data)

        print(f'conversion of plasgraph2 is done {self.output_file}')
        return self.output_file
        