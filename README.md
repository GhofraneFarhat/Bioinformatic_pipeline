# Plaspipe

## Table of Contents
- [Installation](#installation)
- [Pipeline Workflow](#PipelineWorkflow)
- [Pipeline Prerequisites](#PipelinePrerequisites)
- [Usage](#Usage)


Plaspipe is a bioinformatic pipeline for plasmid analysis. This pipeline integrates the bioinformatic tools such as PlASgraph, Platon, Plasbin-flow. It provides a streamlined and efficient workflow for plasmid identification, classification, and binning. The pipeline's flexibility is further enhanced by its support for both [GFA](https://gfa-spec.github.io/GFA-spec/) and [FASTA](https://zhanggroup.org/FASTA/) file formats as input. This versatility allows researchers to analyze both assembled contigs and raw sequencing data, making the pipeline adaptable to various stages of genomic analysis.

## Installation

Plaspipe can be installed from this repository 

```bash
git clone https://github.com/GhofraneFarhat/Bioinformatic_pipeline
```

Plaspipe is composed of this main script:

plaspipe.py is the main script used to process the pipeline.


### Pipeline Workflow

- Input Processing: The pipeline begins by processing the input files (GFA or FASTA).

- Classification: Utilizes classification tools to classify sequences as plasmid or chromosomal.

- Binning: Applies binning tools to group related plasmid sequences into bins.

- Result Integration: Collates and integrates results from all tools into a comprehensive output.

### Pipeline Prerequisites

Plaspipe is written in python (version 3.10.9+).

The Pipeline integrates a suite of cutting-edge bioinformatics tools to provide a comprehensive solution for plasmid analysis. Each tool have it's requires:


- plASgraph requires:
    - NetworkX 2.6.3+
    - Pandas 1.3.5+
    - NumPy 1.21.5+
    - Scikit-learn 0.23.1+
    - Scipy 1.8.0+
    - Biopython 1.79+
    - Matplotlib 3.5.1+
    - TensorFlow 2.8.0+
    - Spektral 1.0.8+

- PlasBin-flow requires:
    - the python module networkX (version 2.7+),
    - the ILP solver gurobi (version 9.1.2+).
    - pandas (version 2.0.0+),
    - matplotlib (version 3.7.0+),
    - seaborn (version 0.12.2+),
    - Biopython (version 1.81+),
    - scipy (version 1.10.1+),
    - numpy (version 1.24.2+).

### Obtaining Gurobi license
PlasBin-flow and PlasBin use the [Gurobi Solver](https://www.gurobi.com/). To use them, a Gurobi license is needed, which is free for academics.

You can obtain and install a free academic license following the instructions provided at [Gurobi: Always Free for Academics](https://www.gurobi.com/academia/academic-program-and-licenses/).

## Running Plaspipe

### Input files

Plaspipe requires the following input files (described in details below) to process the plasmid analysis:

- input file (FASTA/GFA),

- yaml file.

#### input file 

The [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html) file format is used to represent the assembly graph of sequences. It contains information about the contigs and their connections.

The [FASTA](https://zhanggroup.org/FASTA/) file format is used to represent nucleotide sequences. It contains the sequences of the contigs.

#### Yaml file 

To run the Plaspipe, users need to provide a ```YAML configuration file``` that specifies the input data, output directories, and tool parameters. This configuration file allows the pipeline to be flexible and adaptable to various datasets and analysis tools.

```yaml
input:
  path_to_input_gfa: 
  path_to_input_fasta: 

prefix:
outdir_pipeline: 
 
classification:
  name: 
  version: 
  parameters:
  input_format: 
  
  output:
    outdir_classification:

binning:
  name: 
  version:
  parameters: 
  input_format: 
  
  output:
    outdir_binning:

```

## Plaspipe stages

#### Input processing
#### Classification
#### Binning

## Usage

To process an input file, given the input files described above, the command plaspipe.py can be used as follows:

```python
python src/plaspipe.py
       -yf parameter.yaml
       -out output_file
```
where 

- ```parameter.yaml``` is the yaml file of plaspipe configuration
- ```output_file``` is the file name for the output file

