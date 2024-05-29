from Bio import SeqIO
import csv

class FastaToCsv:
    def __init__(self, fasta_file, csv_file):
        self.fasta_file = fasta_file
        self.csv_file = csv_file

    def convert(self):
        with open(self.csv_file, 'w', newline='') as csvfile:
            fieldnames = ['contig_name', 'plasmid_score', 'chromosome_score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()

            for seq_record in SeqIO.parse(self.fasta_file, "fasta"):
                contig_name = seq_record.id

                try:
                    description = seq_record.description.split(' ')[1][0:]
                    print(description)
                except IndexError:
                    print(f"Warning: Unexpected format for contig '{contig_name}'")
                    continue

                plasmid_score = None
                chromosome_score = None

                if description.lower() == '(plasmid)':
                    plasmid_score = 1
                    chromosome_score = 0
                elif description.lower() == '(chromosome)':
                    plasmid_score = 0
                    chromosome_score = 1
                elif description.lower() == '(ambiguous)':
                    plasmid_score = 1
                    chromosome_score = 1
                elif description.lower() == '':
                    plasmid_score = 0
                    chromosome_score = 0
                else:
                    print(f"Warning: Unknown description '{description}' for contig '{contig_name}'")

                writer.writerow({'contig_name': contig_name, 'plasmid_score': plasmid_score, 'chromosome_score': chromosome_score})

        return self.csv_file


"""
conv = FastaToCsv("classify/classification_1716962016.fasta", "resss1.csv")
res = conv.convert()
print(res)
"""