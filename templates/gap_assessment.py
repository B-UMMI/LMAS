#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information on the size of gaps from the filtered assemblies  respective to the reference sequences.
Produces a JSON file with the sample_id, the assembler, and a list of gap sized for all references.

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
import json
from itertools import groupby
import pandas as pd
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "15.12.2020"
__template__ = "GAP_ASSESSMENT-nf"

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

COLUMNS = ['Sample', 'Assembler', 'Reference', 'Reference Length', 'Gap Start', 'Gap End']


def get_gaps(paf_file, ref_name, ref_len):
    """
    Function to process the mapping (*.paf) file for a given reference and output a list with gap sizes from
    the alignment.
    :param paf_file: tabular file with alignment information for an assembler
    :param ref_name: reference name to filter from the paf_filename
    :param ref_len: int with expected reference length
    :return: gaps: list with gap coords in the reference in the assembly for the ref_name reference
    """
    covered_bases_list = []

    with open(paf_file) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                start, end = int(parts[7]), int(parts[8])
                covered_bases_list.append([start, end])

    covered_bases = set()
    for item in sorted(covered_bases_list, key=lambda x: x[0]):
        start, stop = map(int, item[:])
        for base in range(start, stop):
            # Due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth.
            # Therefore, the coordinates are adjusted as follows:
            # [0; ref_len][ref_len+1; 2*ref_len][(2*ref_len)+1; 3*ref_len]
            if base <= ref_len:
                covered_bases.add(int(base))
            elif base <= 2*ref_len:
                covered_bases.add(int(base-ref_len))
            else:
                covered_bases.add(int(base-(2*ref_len)))

    covered_bases = sorted(list(covered_bases))
    gaps = [[s, e]for s, e in zip(covered_bases, covered_bases[1:]) if s+1 < e]  # get list of gap sizes coords
    gap_sizes = [coord[1]-coord[0]-1 for coord in gaps]  # get list of gap sizes
    return gaps, gap_sizes


def main(sample_id, assembler, assembly, mapping, reference):

    all_gap_sizes = []

    df = pd.DataFrame(columns=COLUMNS)

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(open(reference, "r"), lambda line: line[0] == ">"))

    for header in references:
        header_str = header.__next__()[1:].strip().split()[0]
        reference = utils.REFERENCE_DIC[header]
        seq = "".join(s.strip() for s in references.__next__())

        gaps, gap_sizes = get_gaps(mapping, header_str, len(seq) / 3)
        all_gap_sizes.append(gap_sizes)  # for global plot

        # plot gap location per reference per reference
        for coords in gaps:
            df = df.append({'Sample': sample_id, 'Assembler': assembler, 'Reference': reference,
                            'Reference Length': len(seq/3), 'Gap Start': coords[0], 'Gap End': coords[1]},
                           ignore_index=True)

    to_write = {sample_id: {assembler: sorted(all_gap_sizes)}}

    with open("{}_{}_gap_dict.json".format(sample_id, assembler), "w") as fh:
        fh.write(json.dumps(to_write, separators=(",", ":")))

    df.to_csv(sample_id + '_' + assembler + '_gaps.csv')


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE)
