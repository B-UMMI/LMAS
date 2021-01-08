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
from plotly.offline import plot
import plotly.graph_objects as go
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


def check_missassemblies(paf_file):
    """
    :param paf_file: paf files with mapping information
    :return: dictionary with contigs broken with multiple alignment blocks
    """

    missmatch_dict = {"misassembly": {}, "okay": {}}

    with open(paf_file, 'r') as paf_fh:
        for line in paf_fh:
            line = line.split()

            contig, contig_len, query_start, query_end, strand = line[0:5]
            reference, reference_len, target_start, target_end = line[5:9]

            # a non-perfect alignment or a different number of residue matches
            if int(line[11]) != 0 or line[1] != line[9]:

                cigar = line[-1]
                exact_matches, snp, indel = utils.parse_cs(cigar)

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

                if contig in missmatch_dict["misassembly"].keys():
                    missmatch_dict["misassembly"][contig].append(contig_dict)
                else:
                    missmatch_dict["misassembly"][contig] = [contig_dict]
            else:  # contig okay
                contig_dict = {'contig length': contig_len,
                               'query start': int(query_start),
                               'query end': int(query_end),
                               'strand': strand,
                               'reference': reference,
                               'reference length': int(reference_len) / 3,
                               'target start': utils.adjust_reference_coord(int(target_start),
                                                                            int(reference_len) / 3),
                               'target end': utils.adjust_reference_coord(int(target_end),
                                                                          int(reference_len) / 3)
                               }

                if contig in missmatch_dict["okay"].keys():
                    missmatch_dict["okay"][contig].append(contig_dict)
                else:
                    missmatch_dict["okay"][contig] = [contig_dict]

    print(missmatch_dict["misassembly"])

    return missmatch_dict


def calculate_frag_score(assembler, mis_dict):
    """

    :param assembler:
    :param mis_dict:
    :return:
    """

    frag_score_list_mis = []
    frag_score_list_okay = []

    for contig in mis_dict["misassembly"][assembler].keys():
        frag_score = math.log(len(mis_dict["misassembly"][assembler][contig]) / int(mis_dict["misassembly"][assembler][contig][0]['contig length'])) / \
                 math.log(1 / int(mis_dict["misassembly"][assembler][contig][0]['contig length']))
        frag_score_list_mis.append(frag_score)

    for contig in mis_dict["okay"][assembler].keys():
        frag_score = math.log(len(mis_dict["okay"][assembler][contig]) / int(mis_dict["okay"][assembler][contig][0]['contig length'])) / \
                     math.log(1 / int(mis_dict["okay"][assembler][contig][0]['contig length']))
        frag_score_list_okay.append(frag_score)

    print(frag_score_list_mis)
    print(frag_score_list_okay)
    fig = go.Figure()
    fig.add_trace(go.Box(x=frag_score_list_mis,
                         name="misassembled", boxpoints='outliers',
                         boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))
    fig.add_trace(go.Box(x=frag_score_list_okay,
                         name="okay", boxpoints='outliers',
                         boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))

    plot(fig, filename='lala.html', auto_open=True)

    avg_frag_score = sum(frag_score_list_okay + frag_score_list_mis)/len(frag_score_list_okay + frag_score_list_mis)
    print(avg_frag_score)
    return avg_frag_score


def evaluate_misassembled_contigs(mis_dict):
    """

    :param mis_dict:
    :return:
    """
    for assembler in mis_dict.keys():
        missassembled_contigs = 0
        for contig in mis_dict[assembler].keys():
            if len(mis_dict[assembler][contig]) > 1:  # contig broken into multiple alignment blocks
                missassembled_contigs += 1
                #print(contig)
                num_alignment_blocks = len(mis_dict[assembler][contig])
                aligned_bases = 0
                contig_len = int(mis_dict[assembler][contig][0]['contig length'])
                reference = set()
                strands = set()
                blocks_to_order = []
                for alignment_block in mis_dict[assembler][contig]:
                    #print(alignment_block)
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
                #print(order)

                ##### CHECK ORDER OF ASSEMBLY BLOCKS! FROM ORDER; EVALUATE TYPE OF MISASSEMBLY
                #### CHECK STRANDS OF ASSEMBLY BLOCKS
                ### GENETARE DUMMY TEST DATA
        print("missassembled contigs: {}".format(missassembled_contigs))


def main(sample_id, assembler, assembly, mapping):

    # identify contings broken into multiple assembly blocks
    misassembled_contigs = check_missassemblies(mapping)

    #
    avg_frag_score = calculate_frag_score(assembler, misassembled_contigs)
    evaluate_misassembled_contigs(misassembled_contigs["misassembly"])


if __name__ == '__main__':
    #main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING)
    main("mockSample", "Unicycler", "filtered_mockSample_unicycler.fasta", "mockSample_Unicycler.paf")
    #main("mockSample", "GATBMinia", "filtered_mockSample_unicycler.fasta", "mockSample_GATBMiniaPipeline.paf")
    #main("mockSample", "IDBA", "filtered_mockSample_unicycler.fasta", "mockSample_IDBA-UD.paf")
    #main("mockSample", "velvet", "filtered_mockSample_unicycler.fasta", "mockSample_VelvetOptimizer.paf")
