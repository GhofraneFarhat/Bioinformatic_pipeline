import csv
import shutil
import argparse

class Plasgraph2ToCsv:
    def __init__(self, csv_input, csv_output):
        
        self.input_file = csv_input
        self.output_file = csv_output

    def convert(self):
        
        # Open the original CSV file
        with open(self.input_file, 'r') as file:
            reader = csv.reader(file)
            data = list(reader)

        # Swap the columns
        new_data = []
        for row in data[1:]:  # Skip the header row
            new_row = [row[1], row[3],row[2], row[4], row[5]]
            new_data.append(new_row)

        # Write the new CSV file
        with open(self.output_file, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['contig_name','plasmid_score','length','chromosome_score','label'])  # Write the new header
            writer.writerows(new_data)


        return self.output_file