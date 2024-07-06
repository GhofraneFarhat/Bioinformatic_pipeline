import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .plaspipe_utils import process_exception
import sys

# Mandatory fields in GFA contigs and links
GFA_SEQ_KEY = 'Sequence'
GFA_LEN_KEY = 'Length'
GFA_FROM_KEY = 'From'
GFA_FROM_ORIENT_KEY = 'FromOrient'
GFA_TO_KEY = 'To'
GFA_TO_ORIENT_KEY = 'ToOrient'
GFA_OVERLAP_KEY = 'Overlap'

# Conversion of GFA attributes.
# Missing attributes types: B, J
GFA_ATTRIBUTE_TYPE = {
    'i': lambda x: int(float(x)),
    'f': lambda x: float(x),
    'Z': lambda x: str(x),
    'A': lambda x: str(x),
    'H': lambda x: bytes(x)
}

# Open file
def __open_file_read(file_path, gzipped=False):
    """
    Open a file for reading
    """
    if gzipped:
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')
def __open_file_write(file_path, gzipped=False):
    """
    Open a file for writing
    """
    if gzipped:
        return gzip.open(file_path, 'wt')
    else:
        return open(file_path, 'w')





# Conversion from GFA to fasta
def __write_attributes(attributes_dict, keys_to_remove=[], sep=' '):
    """
    Write GFA attributes into a string

    Args:
        - attributes_dict (Dictionary): attribute key -> attribute value
        - keys_to_remove (List(str)): list of keys to not print
        - sep (str): separating string

        Returns:
        (str): list of attributes in format sep.join(key:value)
        attributes with None value are not written
        """
    return sep.join(
        [
            f'{x}:{y}'
            for x, y in attributes_dict.items()
            if x not in keys_to_remove and y is not None
        ]
    )

def __add_attributes(att_data, attributes_list):
    """
    Add attributes to a dictionary

    Args:
        - att_data (List(str)): list of attribute strings in format <key>:<type>:<value>
        - attributes_list (List(str)): list of attribute keys to keep, or ['all']

        Returns:
        - (Dictionary) attribute key -> attribute value
    """
    result = {}
    for att in att_data:
        att_parts = att.split(':')
        if len(att_parts) == 3:
            key, att_type, value = att_parts
            if attributes_list == ['all'] or key in attributes_list:
                try:
                    result[key] = GFA_ATTRIBUTE_TYPE.get(att_type, lambda x: x)(value)
                except Exception as e:
                    print(f"Error processing attribute '{att}': {e}")
        else:
            print(f"Ignoring malformed attribute '{att}'")
    return result

def read_GFA_ctgs(in_file_path, attributes_list, gzipped=False, ctg_fun=lambda x: x, id_fun=lambda x: x):
    """
    Read contigs and their attributes from a GFA file
    """
    result = {}
    with __open_file_read(in_file_path, gzipped) as in_file:
        for gfa_line in [x for x in in_file.readlines() if x[0] == 'S']:
            line = gfa_line.rstrip()
            ctg_data = line.split('\t')
            if len(ctg_data) < 2:
                continue  # Skip lines with fewer than 2 fields
            ctg_id, ctg_seq = ctg_data[1], ctg_data[2]
            ctg_len = len(ctg_seq)
            att_data = [f'{GFA_SEQ_KEY}:Z:{ctg_seq}', f'{GFA_LEN_KEY}:i:{ctg_len}'] + ctg_data[3:]
            result[id_fun(ctg_id)] = ctg_fun(__add_attributes(att_data, attributes_list))
    return result

def write_GFA_to_FASTA(in_GFA_file, out_FASTA_file, in_gzipped, out_gzipped, sep=' '):
    """
    Create a FASTA file from a GFA file
    """
    GFA_ctg_seqs = read_GFA_ctgs(in_GFA_file, attributes_list=['all'], gzipped=in_gzipped)
    ctg_records = [
        SeqRecord(
            Seq(y[GFA_SEQ_KEY]),
            id=x,
            name=x,
            description=f'{x}.GFA {__write_attributes(y, keys_to_remove=[GFA_SEQ_KEY])}'
        )
        for x, y in GFA_ctg_seqs.items()
    ]
    try:
        with __open_file_write(out_FASTA_file, gzipped=out_gzipped) as out_file:
            SeqIO.write(ctg_records, out_file, 'fasta')
    except Exception as e:
        process_exception(f'FASTA/GFA\tWriting {in_GFA_file} to {out_FASTA_file}: {e}')

