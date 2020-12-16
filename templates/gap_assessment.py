#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import os
import math
import re
import json
import pandas as pd
from itertools import groupby
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
    #gap_sizes = [coord[1]-coord[0]-1 for coord in gaps]
    return gaps


def main(sample_id, assembler, assembly, mapping, reference):

    all_gap_distance = []

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(open(reference, "r"), lambda line: line[0] == ">"))

    for header in references:
        header_str = header.__next__()[1:].strip().split()[0]
        seq = "".join(s.strip() for s in references.__next__())

        gaps = get_gaps(mapping, header_str, len(seq) / 3)

        #dist_array = [gap2[0]-gap1[1] for gap2, gap1 in zip(gaps[1:], gaps)]  # this is wrong
        for i in range(1, len(gaps)):
            dist = gaps[i][0] - gaps[i-1][1]
            all_gap_distance.extend(dist)

    to_write = {sample_id: {assembler: sorted(all_gap_distance)}}

    with open("{}_{}_gap_dict.json".format(sample_id, assembler), "w") as fh:
        fh.write(json.dumps(to_write, separators=(",", ":")))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE)
    #main("mockSample", "SKESA", "filtered_mockSample_skesa.fasta", "mockSample_SKESA.paf",
    #     "Zymos_Genomes_triple_chromosomes.fasta")
