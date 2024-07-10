import csv

class PlasBinToCsv:
    def __init__(self, input_file, csv_file):
        self.input_file = input_file
        self.output_file = csv_file

    def convert(self):
        plasmids = {}
        current_plasmid = None

        # Read the input CSV file
        with open(self.input_file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                if row and row[0].startswith('plasmid_'):
                    current_plasmid = row[0].split(';')[0]
                    plasmids[current_plasmid] = {'contigs': []}
                elif row and row[0].startswith('# ') and current_plasmid:
                    parts = row[0].strip('# ').split('\t')
                    if len(parts) >= 3:
                        contig, _, gc_cont = parts[:3]
                        try:
                            plasmids[current_plasmid]['contigs'].append((contig, float(gc_cont)))
                        except ValueError:
                            print(f"Warning: Could not convert GC content to float for contig {contig}")

        print(f"Parsed plasmids: {list(plasmids.keys())}")
    
        if not plasmids:
            print("Error: No plasmids were parsed from the input data.")
            return

        # Write to output CSV file
        try:
            with open(self.output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(['Bin', 'GC_bin', 'Contig'])

                for plasmid, data in plasmids.items():
                    if not data['contigs']:
                        print(f"Warning: No contigs found for {plasmid}")
                        continue

                    gc_values = [gc for _, gc in data['contigs']]
                    min_gc = min(gc_values)
                    max_gc = max(gc_values)
                    gc_bin = f"{min_gc:.2f}-{max_gc:.2f}"
                
                    contig_str = ",".join(f"{contig}:{gc}" for contig, gc in data['contigs'])
                
                    writer.writerow([plasmid, gc_bin, contig_str])

            print(f"CSV file '{self.output_file}' has been created successfully.")
            return self.output_file
        except IOError as e:
            print(f"Error writing to file: {e}")

