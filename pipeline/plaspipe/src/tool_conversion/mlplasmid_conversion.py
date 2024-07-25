import csv



class mlplasmidToCsv:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file


    def convert(self):
    
        header = ['contig_name', 'plasmid_score', 'length', 'chromosome_score', 'label']
    
    
        with open(self.output_file, 'w', newline='') as out_file:
            csvwriter = csv.writer(out_file)
        
        
            csvwriter.writerow(header)
        
        
            with open(self.input_file, 'r') as in_file:
                next(in_file)
            

                for line in in_file:
                    
                    fields = line.strip().split('\t')
            
                    chromosome_score = fields[0].strip('"')
                    plasmid_score = fields[1].strip('"')
                    label = fields[2].strip('"')
                    contig_name = fields[3].split()[1]  
                    length = fields[4]
                
        
                    csvwriter.writerow([contig_name, plasmid_score, length, chromosome_score, label])
        
        return self.output_file