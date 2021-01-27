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
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import pickle
import json
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
    Parses a mapping paf file and stores mapping information into a dictionary for each contig. It separates
    the contigs into two categories (misassembled and okay) depending on the probability of it being misassembled
    (a non-perfect alignment or a different number of residue matches to the contig size).
    :param paf_file: path to the paf files with mapping information
    :return: dictionary with contig mapping info
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

    return missmatch_dict


def get_frag_score(n_blocks, contig_len):
    """
    fragment score calculated by log(n alignment blocks/contig length in bp)/log(1/conting length in bp)
    :param n_blocks: int with number of alignment blocks a contig was divided in
    :param contig_len: int with contign lenght in bp
    :return: float with fragment score
    """
    return math.log(n_blocks/contig_len)/math.log((1/contig_len))


def evaluate_misassembled_contigs(mis_dict):
    """
    Method to evaluate if a contig is misassembed. Recieves a list of possibly misassembled contigs. A contig is
    classified as misassembled if it's broken into multiple alignment blocks. The type of misassembly is classified
    according to the following conditions:
        - Insertion:
        - Deletion:
        - Chimera:
        - Missense:
        - Translocation
    :param mis_dict: dictionary with potentially misassembled contigs
    :return:
    """

    missassembled_contigs = {}
    for contig in mis_dict.keys():
        if len(mis_dict[contig]) > 1:  # contig broken into multiple alignment blocks

            #frag score
            n_blocks = len(mis_dict[contig])
            contig_len = int(mis_dict[contig][0]['contig length'])
            frag_score = (get_frag_score(n_blocks, contig_len))

            # misassembly detection
            reference = set()  # set of references
            strands = set()  # set of strands
            misassembly_list = []
            blocks_coords = []

            blocks_to_order = []  # list with start coordinates of the alignment block in the reference
            for alignment_block in mis_dict[contig]:
                reference.add(alignment_block['reference'])
                strands.add(alignment_block['strand'])
                blocks_to_order.append(alignment_block['query start'])
                blocks_coords.append([alignment_block['query start'], alignment_block['query end']])

            #   -chimera
            if len(reference) > 1:  # if chimera, the rest of the evaluations don't make sense
                misassembly_list.append("chimera {}".format(reference))

            else:
                #   -missense
                if len(strands) > 1 & len(reference) == 1:  # don't count missense if chimera
                    misassembly_list.append("missense")
                #   -multiple alignment blocks to the same reference
                blocks_coords = sorted(blocks_coords, key=lambda x: x[0])
                print(blocks_coords)
                gap_sizes = [blocks_coords[i+1][1] - blocks_coords[i][0] for i in range(0, len(blocks_coords)-1)]
                print(gap_sizes)

                blocks_ordered = sorted(blocks_to_order)
                order = []
                for item in blocks_to_order:
                    order.append(blocks_ordered.index(item))
                print(order)
                if sorted(order) != range(min(order), max(order) + 1):  # check if the order is different from reference
                    # check if lists are inverted:
                    if order == sorted(order, reverse=True):
                        misassembly_list.append("inversion")
                    elif len(order) > 2: # TODO - improve. A might be doing a mistake with over simplefication.. distinguish between inversion, insertion, and translocation
                        #   - translocation
                        misassembly_list.append("translocation")
                #   - insertion
                if all(i > 50 for i in gap_sizes):
                    misassembly_list.append("insertion")
            missassembled_contigs[contig] = {'misassembly': misassembly_list, 'frag_score': frag_score,
                                             'contig length': contig_len, "n blocks": n_blocks}
            print(contig, missassembled_contigs[contig])

    return missassembled_contigs


def main(sample_id, assembler, assembly, mapping):

    # identify contings broken into multiple assembly blocks
    misassembled_contigs = check_missassemblies(mapping)

    mis_contigs = evaluate_misassembled_contigs(misassembled_contigs["misassembly"])
    print(assembler)
    print(mis_contigs)

    report_data = {"sample": sample_id, "assembler": assembler, "misassembled_contigs": len(mis_contigs.keys())}

    x = []
    y = []
    z = []
    for item, value in mis_contigs.items():
        x.append(value['contig length'])
        y.append(value['n blocks'])
        z.append(value['misassembly'])



    df = pd.DataFrame(list(zip(x, y, z)),
                      columns=['Contig Length', 'n blocks', 'Misassembly'])
    # TODO - symbol by df['Misassembly']
    trace = go.Scatter(x=df['Contig Length'], y=df['n blocks'], name=assembler, text=df['Misassembly'],  mode='markers',
                       hovertemplate=
                       "<b>%{text}</b><br><br>" +
                       "Contig Length: %{x:.0f}bp<br>" +
                       "Fragments: %{y:.0}<br>" +
                       "<extra></extra>",)

    with open('{}_{}_trace.pkl'.format(sample_id, assembler), 'wb') as f:
        pickle.dump(trace, f)
    
    with open('{}_{}_contig_lenght.pkl'.format(sample_id, assembler), 'wb') as f:
        pickle.dump(x, f)
    
    with open("{}_{}_misassembly.json".format(sample_id, assembler), "w") as json_report:
        json_report.write(json.dumps(report_data, separators=(",", ":")))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING)
