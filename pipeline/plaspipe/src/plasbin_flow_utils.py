import subprocess
import os
import gzip
import shutil
import sys

from .plaspipe_utils import absolute_path
from .plaspipe_utils import process_arguments
from .plaspipe_utils import gzip_file


def generate_content_gc_file(method_configuration, gfa_file):

    python_executable = sys.executable

    binning_outdir = method_configuration['output']['outdir_binning']
    
    #script_path = method_configuration['plasbin_utils_script']
    plasbin_utils_script = os.path.join(absolute_path(), 'submodules/PlasBin-flow/code/plasbin_utils.py')
            
    out_dir = method_configuration.get('plasbin_out_dir', binning_outdir)
    tmp_dir = method_configuration.get('plasbin_tmp_dir', binning_outdir)
    sample_name = method_configuration['sample_name']
    gc_intervals_file = method_configuration.get('gc_interval_file','')
    if gc_intervals_file == None : 
        gc_intervals_file = ''

    input_file = os.path.join (binning_outdir,'plasbin_flow_sample.txt')

    default_tmp_dir = os.path.join(absolute_path(), 'out_tmp')
    default_out_dir = os.path.join(absolute_path(), 'out/plasbin_flow')
    argument = [sample_name, out_dir, tmp_dir]
    default_arg = ['test1', default_out_dir, default_tmp_dir]
    
    gc_content_argument = process_arguments(argument, default_arg)
    sample_name, out_dir, tmp_dir = gc_content_argument

    #gzip the gfa file for the txt file
    gfa_path = gzip_file(gfa_file, binning_outdir)

    #generate the input file for generate for the gc content
    generate_input_file(input_file, sample_name, gfa_path)

  
    command = [
        python_executable,
        plasbin_utils_script,
        'gc_probabilities', 
        '--input_file', input_file,
        '--out_dir', out_dir,
        '--tmp_dir', tmp_dir,
        '--gc_intervals', gc_intervals_file,
    ]

    print(f'command to run gc_content_file {command}')
    subprocess.run(command, check=True)
    
    #return path_to_gc_content_file

def generate_input_file(output_file, sample = 'ABCD', gfa_path = ''):

    """
    generate un txt file of sample name and GFA file path.

    Args:
        output_file (str): output_file path.
        
    """

    with open(output_file, 'w') as f:

        f.write('sample,gfa\n')
        print(f'we are generating a txt file {output_file}')
        # Écrire le nom d'échantillon et le chemin GFA dans le fichier
        f.write(f'{sample},{gfa_path}\n')


