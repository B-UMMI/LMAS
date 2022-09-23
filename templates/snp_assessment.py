#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information on the location of SNPs from the filtered assemblies  respective to the reference sequences.
Produces a JSON file with the sample_id, the assembler, and a list of SNP locations for all references.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembler``: String with assembler name.
    - e.g.: ``'SPAdes'``
- ``assembly``: fasta file from the assembler (filtered for minimum length size)
    - e.g.: ``'spades.fasta'``
- ``mapping``: paf file of the assembly mapped to the complete triple reference genomes
    - e.g.: ``'spades.paf' ``
- ``reference``: paf file to the complete triple reference genomes
    - e.g.: ``'references.fasta' ``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import re
from itertools import groupby
from numpy.lib.function_base import append
import pandas as pd
import json
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "16.02.2021"
__template__ = "SNP_ASSESSMENT-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    MAPPING = '$mapping'
    REFERENCE = '$reference'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("MAPPING: {}".format(MAPPING))
    logger.debug("REFERENCE: {}".format(REFERENCE))

COLUMNS = ['Sample', 'Assembler', 'Reference', 'Reference Length', 'SNP Location']


def get_position(start, end, cigar):
    """
    """
    #parse cigar components - exact matches, insertions, deletions and substitutions
    matches = re.findall(r'([:=*+-])(\\d+|[A-Za-z]+)', cigar)
    
    coord = start
    substitution = []
    for match in matches:
        if match[0] == ':':
            if match[1] == 'Z':
                continue
            if int(match[1]): #exact matches to the ref
                coord += int(match[1])
        if match[0] == '-': # if deletion, coords in the ref need to be ajusted
            coord += len(match[1])
        #if match[0] == '+': # insertions to the reference don't matter
            #coord += len(match[1])
        if match[0] == '*': 
            substitution = match[1]
            coord += 1
            yield coord, substitution


def get_snps(paf_file, ref_name, ref_len, sample_id, assembler):
    """
    Function to process the mapping (*.paf) file for a given reference and output a list with snp locations from
    the alignment.
    :param paf_file: tabular file with alignment information for an assembler
    :param ref_name: reference name to filter from the paf_filename
    :param ref_len: int with expected reference length
    :return: snps: list with snp coords in the reference in the assembly for the ref_name reference
    """
    snps = []
    tsv_report = open("{}_{}_{}_substitutions.tsv".format(sample_id, assembler, ref_name), "w")

    with open(paf_file) as paf:
        for line in paf:
            parts = line.strip().split()
            if parts[5] == ref_name:
                if parts[4] == '+':
                    start, end = int(parts[7]), int(parts[8])
                else:
                    start, end = int(parts[8]), int(parts[7])
                cigar = parts[-1]
                if len(re.findall(r'\\*', cigar)) > 0:
                    snps_iterator = get_position(start, end, cigar)
                    for snp in snps_iterator:
                        tsv_report.write('\\t'.join([str(utils.adjust_reference_coord(snp[0], ref_len)), str(snp[1][0]), str(snp[1][1])]) + '\\n')
                        snps.append((utils.adjust_reference_coord(snp[0], ref_len), snp[1]))
                else:
                    continue

    tsv_report.close()
    return snps


def main(sample_id, assembler, assembly, mapping, reference):

    df = pd.DataFrame(columns=COLUMNS)

    reference_report = {"sample": sample_id,
                    "assembler": assembler, "reference": {}}

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(open(reference, "r"), lambda line: line[0] == ">"))

    for header in references:
        reference_name = header.__next__()[1:].strip()
        print(reference_name)
        seq = "".join(s.strip() for s in references.__next__())
        snps = get_snps(mapping, reference_name, len(seq) / 3, sample_id, assembler)

        # plot gap location per reference per reference
        for snip_info in snps:
            coord = snip_info[0]
            substitution = '{}->{}'.format(snip_info[1][0], snip_info[1][1])
            #print(coord, substitution)
            df = df.append({'Sample': sample_id, 'Assembler': assembler, 'Reference': reference_name,
                            'Reference Length': len(seq)/3, 'SNP Location': coord, 'Substitution Type': substitution}, ignore_index=True)
            
            if reference_name not in reference_report['reference'].keys():
                reference_report['reference'][reference_name] = {"snps": 1}
            else:
                reference_report['reference'][reference_name]["snps"] += 1

    # Write files for report
    with open("{}_{}_snps.json".format(sample_id, assembler), "w") as json_report:
        json_report.write(json.dumps(reference_report, separators=(",", ":")))
    
    df.to_csv(sample_id + '_' + assembler + '_snps.csv')


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE)
