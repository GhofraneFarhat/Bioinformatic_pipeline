import csv
import shutil
import argparse

class CsvToCsv:
    def __init__(self, csv_input, csv_output):
        
        self.input_file = csv_input
        self.output_file = csv_output

    def convert(self):
        
        shutil.copy(self.input_file, self.output_file)
        return self.output_file


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='CSV to CSV')
    parser.add_argument('csv_file_input', type=str, help='Path to the input csv file')
    parser.add_argument('csv_file_output', type=str, help='Path to the output csv file')

    args = parser.parse_args()

    conversion = CsvToCsv(args.csv_file_input, args.csv_file_output)
    res = conversion.convert()