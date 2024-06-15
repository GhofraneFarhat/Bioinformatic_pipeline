import subprocess

def generate_content_gc_file(method_configuration, gfa_path):

    binning_outdir = method_configuration['outdir_binning']
    
    script_path = method_configuration['plasbin_utils_script']
    out_dir = method_configuration.get('plasbin_out_dir', binning_outdir)
    tmp_dir = method_configuration.get('plasbin_tmp_dir', binning_outdir)
    sample_name = method_configuration['sample_name']

    input_file = os.path.join (method_configuration['outdir_binning'],'plasbin_flow_sample.txt')
    
    generate_input_file(input_file, sample_name, gfa_path)

    command = ["python", script_path, "gc_probabilities", "--inpt_file", input_file, "--out_dir", out_dir, "--tmp_dir", tmp_dir]

    subprocess.run(command, check=True)

def generate_input_file(gfa_path, sample = 'ABCD'):

    """
    Génère un fichier texte contenant des noms d'échantillon et des chemins GFA.

    Args:
        output_file (str): Chemin du fichier de sortie.
        num_samples (int, optional): Nombre d'échantillons à générer. Défaut à 10.
    """

    with open(output_file, 'w') as f:
        for _ in range(num_samples):
            

            # Écrire le nom d'échantillon et le chemin GFA dans le fichier
            f.write(f'{sample},{gfa_path}\n')