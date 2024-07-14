import csv


class PlasForestToCsv:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file

    def convert(self):
        with open(self.input_file, 'r') as in_file, open(self.output_file, 'w', newline='') as out_file:
            
            csv_writer = csv.writer(out_file)
        
            # Write the header
            csv_writer.writerow(['contig_name', 'plasmid_score', 'chromosome_score'])
        
            
            next(in_file) #skip first line
        
            
            for line in in_file:
                id, prediction = line.strip().split(',')
            
                if prediction.lower() == 'plasmid':
                    plasmid_score = 1
                    chromosome_score = 0

                elif prediction.lower() == 'ambiguous':
                    plasmid_score = 1
                    chromosome_score = 1

                elif prediction.lower() == '':
                    plasmid_score = 0
                    chromosome_score = 0

                else:  
                    plasmid_score = 0
                    chromosome_score = 1
            
                # Write the transformed row
                csv_writer.writerow([id, plasmid_score, chromosome_score])

        print(f"CSV file '{self.output_file}' has been created successfully.")
        return self.output_file



