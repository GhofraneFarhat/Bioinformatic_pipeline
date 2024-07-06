import sys #handle command-line arguments
import re
from collections import defaultdict

def classify(input_fasta, output_fasta):
    classifications = defaultdict(str)

    with open(input_fasta, 'r') as file:
        for line in file:
            if line.startswith('>'):
                contig_name = line.strip('>\n')

                # with names
                if 'chromosome' in contig_name.lower():
                    classifications[contig_name] = 'chromosome'
                elif 'plasmid' in contig_name.lower():
                    classifications[contig_name] = 'plasmid'
                else:
                    classifications[contig_name] = 'ambiguous'

    with open(output_fasta, 'w') as output_file:
        for contig, classification in classifications.items():
            output_file.write(f'>{contig} ({classification})\n')

            # result
            with open(input_fasta, 'r') as input_file:
                found_contig = False
                for seq_line in input_file:
                    if seq_line.startswith(f'>{contig}'):
                        found_contig = True
                    elif found_contig and not seq_line.startswith('>'):
                        output_file.write(seq_line)
                    elif found_contig and seq_line.startswith('>'):
                        break

    return output_fasta

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    classify(input_fasta, output_fasta)
    print(f"Results written to {output_fasta}")