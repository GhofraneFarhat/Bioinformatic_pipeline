import csv
import os

def csv_to_tsv (csv_file, path_to_csv, path_to_tsv):
   
    tool_result = os.path.join(path_to_tsv, 'plasbin-flow_output.tsv')
    csv.writer(open(tool_result, 'w+'), delimiter='\t').writerows(csv.reader(open(path_to_csv)))
    print(tool_result)
