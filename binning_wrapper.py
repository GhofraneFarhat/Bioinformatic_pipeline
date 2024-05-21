import os
import subprocess
from Bio import SeqIO


class BinningWrapper:


#initialisation
    def __init__(self, input_file, output_file, method_config, repo_path, pipeline_data):
        self.input_path = input_file
        self.output_path = output_file
        self.method_configuration = method_config
        self.repo_path = repo_path
        #self.bin_dir = temp_dir
        self.pipeline_data = pipeline_data
        #the path to save the output of the tool        


#to run the tool frol the file load.py
    def run(self):

        bin_format = self.method_configuration['input_format']  # bin_format for the format of the input of the binning tool
        #input_bin = self.pipeline_data.conversion(bin_format) #input_bin: the path to input file

        # getting the tool parameter using the yamk file
        tool_parameters = self.method_configuration['parameters'] # the parameter of the binning tool

        # Run the binning tool : need to see more about the parameters
        self.run_binning()

        # Update the PipelineData instance with the results
        self.update_pipeline_data()


#the function to run the binning tool
    def run_binning(self):
    # Extract the tool_command from the method_config dictionary
        tool_command = self.method_configuration.get("command")
        #conv = Conversion(self.input_path)
        #input_file = conv.conversion('fasta')
    # Split the tool_command into a list
        command_parts = tool_command.split()

    # Construct the full command with input_file_path and output_file
        full_command = command_parts + [(self.input_path), self.output_path]

    # Change the current working directory to the repository path
        os.chdir(self.repo_path)

    # Run the command using subprocess
        subprocess.run(full_command, shell=False)


#the function to update the pipeline instance

    # update the pipeline data using the binning tool
    def update_pipeline_data(self):
        print("let's update")
        #output_bin = os.path.join(self.bin_dir, 'bin.fasta')
        output_bin = self.output_path

        with open(output_bin, 'r') as binning_result:
            #output csv
            if self.method_configuration['output_format'] == 'csv':
                self.update_pipeline_data_from_csv(binning_result)
            #output tsv
            elif self.method_configuration['output_format'] == 'tsv':
                self.update_pipeline_data_from_tsv(binning_result)
            #output txt
            elif self.method_configuration['output_format'] == 'txt':
                self.update_pipeline_data_from_txt(binning_result)
            #output fasta
            elif self.method_configuration['output_format'] == 'fasta':
                self.update_pipeline_data_from_fasta(binning_result)
            else:
                print("Unsupported output format")

                #we can add other tools or find a solution for a general format of the output

# update from a csv file
    def update_pipeline_data_from_csv(self, binning_result):
        # Skip the first line
        next(binning_result_result)
        
        for line in binning_result:
            contig_name, bin_id = line.strip().split(',')

            contig_list = self.pipeline_data.get_bin(bin_id)
            if contig_list is None:
                self.pipeline_data.set_bins(bin_id, [contig_name])
            else:
                contig_list.append(contig_name)
                self.pipeline_data.set_bins(bin_id, contig_list)


    #we need to handle if in the conversation of the file we lose a contig cause of an error 


 # update from a tsv file 
    def update_pipeline_data_from_tsv(self, binning_result):
        
        for line in binning_result:
            contig_name, bin_id = line.strip().split('\t')

            contig_list = self.pipeline_data.get_bin(bin_id)
            if contig_list is None:
                self.pipeline_data.set_bins(bin_id, [contig_name])
            else:
                contig_list.append(contig_name)
                self.pipeline_data.set_bins(bin_id, contig_list)

# update from a txt file
    def update_pipeline_data_from_txt(self, binning_result):
        
        pass

#update from fasta file
    def update_pipeline_data_from_fasta(self, binning_result):
        
        for seq_record in SeqIO.parse(binning_result, "fasta"):
            contig_name = seq_record.id
            bin_id = seq_record.description.split()
            #print(description[1])
            
            contig_list = self.pipeline_data.get_bin(bin_id)
            if contig_list is None:
                self.pipeline_data.set_bins(bin_id, [contig_name])
            else:
                contig_list.append(contig_name)
                self.pipeline_data.set_bins(bin_id, contig_list)            

        #we need to handle if in the conversation of the file we lose a contig cause of an error 