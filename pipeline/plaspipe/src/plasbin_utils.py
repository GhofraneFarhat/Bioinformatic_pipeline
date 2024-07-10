import subprocess
import os
import gzip
import shutil
import glob
import pandas as pd
import logging
import sys 

from .plaspipe_utils import absolute_path
from .plaspipe_utils import clean_files
from .plaspipe_utils import log_file
from .plaspipe_utils import run_cmd


def generate_seed_contigs(gfa_file, mapping_file, method_configuration):

    binning_outdir = method_configuration['output']['outdir_binning']
    out_dir = method_configuration.get('plasbin_out_dir', binning_outdir)

    
    seed_contigs_script = os.path.join(absolute_path(), 'submodules/PlasBin/code/generate_seeds.py')
       

    command = ["python", seed_contigs_script, "--ag", gfa_file, "--map", mapping_file, "--out", out_dir]


    print(f'command to run seed_contigs_file {command}')
    subprocess.run(command, check=True)

    #return path_to_seed_contigs

def run_blast6(query_file, db_file, mappings_file):
    db_prefix = f'{db_file}.db'
    logging.info(f'ACTION\tcompute blast database for {db_file}: {db_prefix}')
    
    
    cmd_makeblastdb = [
        'makeblastdb',
        '-in', db_file,
        '-dbtype', 'nucl',
        '-out', db_prefix
    ]
    try:
        _ = run_cmd(cmd_makeblastdb)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running makeblastdb: {e}")
        raise

    logging.info(f'ACTION\tmap {query_file} to {db_prefix}')        
    cmd_megablast = [
        'blastn', '-task', 'megablast',
        '-query', query_file,
        '-db', db_prefix,
        '-out', mappings_file,
        '-outfmt', '6'
    ]
    try:
        _ = run_cmd(cmd_megablast)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error running blastn: {e}")
        raise

    print('BLAST completed successfully')
    log_file(mappings_file)
    db_files = glob.glob(f'{db_prefix}.n*')
    clean_files(db_files)




