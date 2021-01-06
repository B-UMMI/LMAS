#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
For each contig, this script will evaluate if it's a missassembly and classy it into these main types:
    * Insertion
    * Deletion
    * Missense
    * Inversion
    * Chimera
    * Translocation
    * ...

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the metagenomic assembly files (ending in *.fasta)
  * Path to the mapped contigs to the triple reference genomes (ending in *.paf)

The triple bacterial reference files for the zymos mock community are available at
"../../data/references/Zymos_Genomes_triple_chromosomes.fasta"

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import math
import os
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "06.01.2021"
__template__ = "MISASSEMBLY-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    MAPPING = '$mapping'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("MAPPING: {}".format(MAPPING))

def check_missassemblies(assembler, paf_file):
    """
    :param assembler: string with assembler name
    :param paf_file: paf files with mapping information
    :return: dictionary with contigs broken with multiple alignment blocks
    """

    missmatch_dict = {}

    with open(paf_file, 'r') as paf_fh:
        for line in paf_fh:
            line = line.split()

            contig, contig_len, query_start, query_end, strand = line[0:5]
            reference, reference_len, target_start, target_end = line[5:9]

            # a non-perfect alignment or a different number of residue matches
            if int(line[11]) != 0 or line[1] != line[9]:

                cigar = line[-1]
                exact_matches, snp, indel = utils.parse_cs(cigar)  # TODO - gap size to be adjusted by param

                contig_dict = {'contig length': contig_len,
                               'query start': int(query_start),
                               'query end': int(query_end),
                               'strand': strand,
                               'reference': reference,
                               'reference length': int(reference_len)/3,
                               'target start': utils.adjust_reference_coord(int(target_start),
                                                                            int(reference_len)/3),
                               'target end': utils.adjust_reference_coord(int(target_end),
                                                                          int(reference_len)/3),
                               'exact matches': exact_matches,
                               'snp': snp,
                               'indels': indel}

                if assembler not in missmatch_dict.keys():
                    missmatch_dict[assembler] = {contig: [contig_dict]}
                else:
                    if contig in missmatch_dict[assembler].keys():
                        missmatch_dict[assembler][contig].append(contig_dict)
                    else:
                        missmatch_dict[assembler][contig] = [contig_dict]

    return missmatch_dict


def evaluate_misassembled_contigs(mis_dict):
    """

    :param mis_dict:
    :return:
    """
    for assembler in mis_dict.keys():
        for contig in mis_dict[assembler].keys():
            if len(mis_dict[assembler][contig]) > 1:  # contig broken into multiple alignment blocks
                print(contig)
                num_alignment_blocks = len(mis_dict[assembler][contig])
                aligned_bases = 0
                contig_len = int(mis_dict[assembler][contig][0]['contig length'])
                reference = set()
                strands = set()
                blocks_to_order = []
                for alignment_block in mis_dict[assembler][contig]:
                    print(alignment_block)
                    aligned_bases += (alignment_block['query end'] - alignment_block['query start'])
                    reference.add(alignment_block['reference'])
                    strands.add(alignment_block['strand'])
                    blocks_to_order.append(alignment_block['query start'])
                    # get order
                if not math.isclose(aligned_bases, contig_len, rel_tol=50):  # TODO - Hardcoded!
                    print("has gaps")
                if len(reference) > 1:
                    print("Different references! {}".format(reference))
                if len(strands) > 1:
                    print("Different strands! {}".format(strands))
                # get order of blocks
                blocks_ordered = sorted(blocks_to_order)
                order = []
                for item in blocks_to_order:
                    order.append(blocks_ordered.index(item))
                print(order)

                ##### CHECK ORDER OF ASSEMBLY BLOCKS! FROM ORDER; EVALUATE TYPE OF MISASSEMBLY
                #### CHECK STRANDS OF ASSEMBLY BLOCKS
                ### GENETARE DUMMY TEST DATA


def main(sample_id, assembler, assembly, mapping):

    # identify contings broken into multiple assembly blocks
    misassembled_contigs = check_missassemblies(mapping)
    evaluate_misassembled_contigs(misassembled_contigs)


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING)
