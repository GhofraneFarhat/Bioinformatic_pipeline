import csv
from collections import defaultdict
import sys

def bin(input_file, output_file):

    contig_bins = defaultdict(str)
    
    with open(input_file, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
                
                
            if len(fields) >= 3 and fields[0] == 'S':
                contig_name = fields[1]
                contig_sequence = fields[2]
                try:
                    contig_length = len(contig_sequence)
                except ValueError:
                    print(f"Warning: Invalid contig sequence '{contig_sequence}' encountered. Skipping this line.")
                    continue
                
                # length
                if contig_length < 1000:
                    bin_name = 'Bin 1'
                elif contig_length < 5000:
                    bin_name = 'Bin 2'
                elif contig_length < 10000:
                    bin_name = 'Bin 3'
                else:
                    bin_name = 'Bin 4'
                
                contig_bins[contig_name] = bin_name
    
    # csv file
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Contig', 'Bin'])
        for contig, bin_name in contig_bins.items():
            writer.writerow([contig, bin_name])

    return output_file


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_gfa> <output_csv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    bin(input_file, output_file)
    print(f"Results written to {output_file}")