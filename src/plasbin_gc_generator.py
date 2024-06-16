import subprocess
import os
import gzip
import shutil


def generate_content_gc_file(method_configuration, gfa_file):

    binning_outdir = method_configuration['output']['outdir_binning']
    
    script_path = method_configuration['plasbin_utils_script']
    out_dir = method_configuration.get('plasbin_out_dir', binning_outdir)
    tmp_dir = method_configuration.get('plasbin_tmp_dir', binning_outdir)
    sample_name = method_configuration['sample_name']
    gc_intervals_file = method_configuration.get('gc_interval_file','')
    if gc_intervals_file == None : 
        gc_intervals_file = ''

    input_file = os.path.join (binning_outdir,'plasbin_flow_sample.txt')

    #gzip the gfa file for the txt file
    gfa_path = gzip_file(gfa_file)

    #generate the input file for generate for the gc content
    generate_input_file(input_file, sample_name, gfa_path)

    command = ["python", script_path, "gc_probabilities", "--input_file", input_file, "--out_dir", out_dir, "--tmp_dir", tmp_dir, "--gc_intervals", gc_intervals_file]

    print(f'command to run gc_content_file {command}')
    subprocess.run(command, check=True)

    #return path_to_gc_content_file

def generate_input_file(output_file, sample = 'ABCD', gfa_path = ''):

    """
    Génère un fichier texte contenant des noms d'échantillon et des chemins GFA.

    Args:
        output_file (str): Chemin du fichier de sortie.
        num_samples (int, optional): Nombre d'échantillons à générer. Défaut à 10.
    """

    with open(output_file, 'w') as f:

        f.write('sample,gfa\n')
        print(f'we are generating a txt file {output_file}')
        # Écrire le nom d'échantillon et le chemin GFA dans le fichier
        f.write(f'{sample},{gfa_path}\n')

def gzip_file(file_path):
    """
    Gzips the given file and creates a new file with '.gz' extension.
    
    Args:
        file_path (str): Path to the file to be gzipped.
    """
    # Get the file name and directory
    directory, file_name = os.path.split(file_path)
    
    # Create the gzipped file path
    gzipped_file_path = os.path.join(directory, file_name + '.gz')
    
    # Open the input file in binary mode
    with open(file_path, 'rb') as input_file:
        # Create the gzipped file and write the compressed data
        with gzip.open(gzipped_file_path, 'wb') as gzipped_file:
            shutil.copyfileobj(input_file, gzipped_file)

    return gzipped_file_path
    
    print(f"File '{file_path}' has been gzipped to '{gzipped_file_path}'")