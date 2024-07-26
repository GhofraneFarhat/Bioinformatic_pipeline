#!/usr/bin/env python3

"""
Command Generator for Bioinformatics Tools

This module provides functionality to generate command-line instructions for various
bioinformatics tools used in plasmid analysis and classification.

"""

# Import statements
from .plasbin_flow_utils import generate_content_gc_file, gzip_file
from .plaspipe_utils import (
    csv_to_tsv, 
    check_file, 
    process_exception, 
    process_error, 
    process_arguments,
    absolute_path, 
    conversion_gfa_fasta, 
    copy_file, 
    csv_to_tab, 
    update_contig_names
)
import os
import logging
import subprocess
import sys

def get_command(method_config, input_file, output_file, plasmid_scores_file=""):
    """
    Generate the command for running various bioinformatics tools based on the method configuration.

    Args:
        method_config (dict): Configuration dictionary for the method.
        input_file (str): Path to the input file.
        output_file (str): Path to the output file.
        plasmid_scores_file (str): Path to the plasmid scores file (optional).

    Returns:
        list or str: The generated command to run the tool.
    """
    try:
        # Extract tool information from method_config
        tool_name = method_config['name']
        tool_version = method_config['version']
        
        logging.info(f"Generating command for {tool_name} version {tool_version}")
        logging.debug(f"Input file: {input_file}")
        logging.debug(f"Output file: {output_file}")

        python_executable = sys.executable

        # Generate command based on tool and version
        if tool_name == "plASgraph2" and tool_version == "1.0.0":
            return generate_plasgraph2_command(method_config, input_file, output_file)
        elif tool_name == "platon" and tool_version == "1.0.0":
            return f'platon {input_file} {output_file}'
        elif tool_name == "plasbin_flow" and tool_version == "1.0.2.2":
            return generate_plasbin_flow_command(method_config, input_file, output_file, plasmid_scores_file)
        elif tool_name == "PlasForest" and tool_version == "1.4.0":
            return generate_plasforest_command(input_file, output_file)
        elif tool_name == "RFPlasmid" and tool_version == "1.0.0":
            return generate_rfplasmid_command(input_file, output_file)
        elif tool_name == "gplas2" and tool_version == "1.1.0":
            return generate_gplas2_command(method_config, input_file, output_file, plasmid_scores_file)
        elif tool_name == "mlplasmid" and tool_version == "2.2.0":
            return generate_mlplasmid_command(method_config, input_file, output_file)
        else:
            raise ValueError(f"unsupported tool: {tool_name} or version: {tool_version}")

    except KeyError as e:
        process_error(f"missing key in method_config: {e}")
    except Exception as e:
        process_exception(f"unexpected error in get_command: {e}")

def generate_plasgraph2_command(method_config, input_file, output_file):
    """Generate command for plASgraph2"""
    model = method_config.get('model', os.path.join(absolute_path(), 'submodules/plASgraph2/model/ESKAPEE_model'))
    plasgraph2_script = os.path.join(absolute_path(), 'submodules/plASgraph2/src/plASgraph2_classify.py')
    return [
        sys.executable,
        plasgraph2_script,
        'gfa',
        input_file,
        model,
        output_file
    ]

def generate_plasbin_flow_command(method_config, input_file, output_file, plasmid_scores_file):
    """Generate command for plasbin_flow"""
    path_to_outdir, output_file = os.path.split(output_file)
    sample_name = method_config['sample_name']
    log_file = method_config.get('log_file', 'plasbin_flow.log')
    assembler = method_config['assembler']
    seed_score = method_config['seed_score']
    seed_len = method_config['seed_len']
    gc_intervals = method_config['gc_interval_file']
    out_gc_file = method_config['plasbin_out_dir']

    # Set default values
    default_gc_intervals = os.path.join(absolute_path(), 'submodules/PlasBin-flow/example/default/gc_intervals.txt')
    default_out_dir = os.path.join(absolute_path(), 'out/plasbin_flow')
    #script for plasbin_flow
    plasbin_flow_script = os.path.join(absolute_path(), 'submodules/PlasBin-flow/code/plasbin_flow.py')

    # Process arguments
    plasbin_flow_argument = process_arguments(
        [sample_name, assembler, seed_score, seed_len, gc_intervals, out_gc_file],
        ['test1', 'unicycler', 0.58, 2650, default_gc_intervals, default_out_dir]
    )
    sample_name, assembler, seed_score, seed_len, gc_intervals, out_gc_file = plasbin_flow_argument

    # Generate GC content file
    try:
        generate_content_gc_file(method_config, input_file)
        gc_content_file = os.path.join(out_gc_file, f"{sample_name}.gc.tsv")
        check_file(gc_content_file)
        logging.info(f"GC content file is saved to {gc_content_file}")
    except Exception as e:
        process_error(f"Error generating or accessing GC content file: {e}")

    # Gzip input file
    try:
        input_file = gzip_file(input_file)
        check_file(input_file)
        logging.info(f"Gzipped file {input_file}")
    except Exception as e:
        process_error(f"Error gzipping input file: {e}")

    # Convert CSV plasmid score file to TSV
    try:
        classification_score_result = csv_to_tsv(plasmid_scores_file)
        check_file(classification_score_result)
    except Exception as e:
        process_error(f"Error converting plasmid scores file to TSV: {e}")

    return [
        sys.executable,
        plasbin_flow_script,
        '-ag', input_file,
        '-gc', gc_content_file,
        '-score', classification_score_result,
        '-out_dir', path_to_outdir,
        '-out_file', output_file,
        '-log_file', log_file,
        '-seed_len', str(seed_len),
        '-seed_score', str(seed_score),
        '-assembler', assembler,
        '-gc_intervals', gc_intervals
    ]

def generate_plasforest_command(input_file, output_file):
    """Generate command for PlasForest"""
    bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/PlasForest.sh')
    os.chmod(bash_script, 0o755)
    return [bash_script, input_file, output_file]

def generate_rfplasmid_command(input_file, output_file):
    """Generate command for RFPlasmid."""
    project_path = absolute_path()
    path_to_outdir, _ = os.path.split(output_file)
    bash_script = os.path.join(project_path, 'pipeline/plaspipe/src/RFPlasmid.sh')
    #create the fasta folder for rfplasmid input
    fasta_folder = os.path.join(project_path, "fasta_folder")
    os.makedirs(fasta_folder, exist_ok=True)
    copy_file(input_file, fasta_folder)

    os.chmod(bash_script, 0o755)
    return [bash_script, fasta_folder, path_to_outdir, project_path]

def generate_gplas2_command(method_config, input_file, output_file, plasmid_scores_file):
    """Generate command for gplas2"""
    bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/gplas2.sh')
    min_length, output_name = process_arguments(
        [method_config['min_length'], method_config['output_name']],
        [100, 'gplas']
    )
    os.chmod(bash_script, 0o755)

    try:
        classification_result = csv_to_tab(plasmid_scores_file)
        check_file(classification_result)
        classification_result_file = update_contig_names(classification_result, input_file)
        check_file(classification_result_file)
    except Exception as e:
        process_error(f"Error processing files for gplas2: {e}")

    return [
        "bash",
        bash_script,
        input_file,
        classification_result_file,
        output_name,
        str(min_length)
    ]

def generate_mlplasmid_command(method_config, input_file, output_file):
    """Generate command for mlplasmid."""
    path_to_outdir, _ = os.path.split(output_file)
    threshold = method_config.get('threshold', 0.7)
    species = method_config['species']

    mlplasmid_species = ['Enterococcus faecium', 'Klebsiella pneumoniae', 'Escherichia coli']
    if species not in mlplasmid_species:
        process_error(f"Unsupported species '{species}'. Supported species are: {', '.join(mlplasmid_species)}")

    input_file = gzip_file(input_file, path_to_outdir)
    bash_script = os.path.join(absolute_path(), 'pipeline/plaspipe/src/mlplasmid.sh')
    os.chmod(bash_script, 0o755)

    return [bash_script, input_file, output_file, str(threshold), species]

