import subprocess
import os
import gzip
import shutil

from .plaspipe_utils import absolute_path


def generate_seed_contigs(gfa_file, mapping_file, method_configuration):

    binning_outdir = method_configuration['output']['outdir_binning']
    out_dir = method_configuration.get('plasbin_out_dir', binning_outdir)

    
    seed_contigs_script = os.path.join(absolute_path(), 'submodules\PlasBin\code\generate_seeds.py')
       

    command = ["python", seed_contigs_script, "--ag", gfa_file, "--map", mapping_file, "--out", out_dir]


    print(f'command to run seed_contigs_file {command}')
    subprocess.run(command, check=True)

    #return path_to_seed_contigs