#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
For each contig, this script will evaluate if it's a missassembly and classy it into these main types:
    * Insertion
    * Deletion
    * Inversion
    * Chimera
    * Translocation
    * Duplication

Expected input
--------------
This script takes the following arguments (in this order):
  * Path to the metagenomic assembly files (ending in *.fasta)
  * Path to the mapped contigs to the triple reference genomes (ending in *.paf)

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import math
import os
from re import S
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import pickle
import json
from plotly.validators.scatter.marker import SymbolValidator
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.2"
__build__ = "01.06.2021"
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


def parse_paf(paf_file):
    """
    Parses a mapping paf file and stores mapping information into a dictionary for each contig.
    :param paf_file: path to the paf files with mapping information.

    :return: dictionary with contig mapping info
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

            if contig in missmatch_dict.keys():
                missmatch_dict[contig].append(contig_dict)
            else:
                missmatch_dict[contig] = [contig_dict]

    return missmatch_dict


def filter_dict(paf_dict):
    """
    Filters contig dictionary to contain contigs that were broken into more
    than one alignment block.

    Parameters
    ----------
    paf_dict : dict
        dictionary with contig mapping info

    Returns
    -------
    dict
        filtered dictionary containing misassembled contigs
    """

    filtered_dict = {}
    for contig in paf_dict.keys():
        if len(paf_dict[contig]) > 1:
            filtered_dict[contig] = paf_dict[contig]

    return filtered_dict


def classify_misassembled_contigs(mis_dict):
    """
    Method to classify a misassembled contig. Recieves a dictionary of misassembled contigs.
    A contig is classified as misassembled if it's broken into multiple alignment blocks.
    The misassembly is classified into to the following categories:
        - Chimera: blocks align to more than one reference
        - Inversion: blocks align to more than one strand for the same reference
        - Translocation: Blocks align to different regions in the same reference
        - Inconsistency: Catch all?
    :param mis_dict: filtered dictionary containing misassembled contigs
    :return: dictionary classified misassembled contigs
    """

    missassembled_contigs = {}
    for contig in mis_dict:

        # list with misassembly classification (in case of multiple events)
        misassembly_list = []
        n_blocks = len(mis_dict[contig])  # number of alignment blocks
        # total contig lenght
        contig_len = int(mis_dict[contig][0]['contig length'])

        # ---- misassembly detection datastructures ---- #
        reference = set()  # set of references
        strands = set()  # set of strands
        # list of tuples with alignment block coordenates in contig
        blocks_coords_in_contig = []
        # list of tuples with alignment block coordenates in contig
        blocks_coords_in_reference = []
        ref_length = set()
        distances_between_blocks_ref = []
        distance_between_blocks_contig = []

        # ---- fill out datastructures ---- #
        for alignment_block in mis_dict[contig]:
            reference.add(alignment_block['reference'])
            strands.add(alignment_block['strand'])
            blocks_coords_in_contig.append(
                [alignment_block['query start'], alignment_block['query end']])
            blocks_coords_in_reference.append(
                [alignment_block['target start'], alignment_block['target end']])
            ref_length.add(alignment_block['reference length'])

        # sort blocks order in reference and in contigs on start position
        sorted_blocks_coords_in_reference = sorted(
                blocks_coords_in_reference, key=lambda x: x[0])
        sorted_blocks_coords_in_contig = sorted(
                blocks_coords_in_contig, key=lambda x: x[0])
        
        # check blocks order
        order_in_ref = []
        for item in blocks_coords_in_reference:
            order_in_ref.append(sorted_blocks_coords_in_reference.index(item))

        #### MISASSEMBLY CLASSIFICATION ALGORYTHM ####
        #   A. chimera
        if len(reference) > 1:  # if chimera, the rest of the evaluations don't make sense
            misassembly_list.append("chimera {}".format(reference))

        #   B. only one reference
        else:

            # B1. inversion
            if len(strands) > 1:
                misassembly_list.append("inversion")

            # B2. Check distance between alignment blocks

            # B2.1 In Reference
            for i in range(0, len(sorted_blocks_coords_in_reference)-1):
                distances_between_blocks_ref.append(
                    sorted_blocks_coords_in_reference[i+1][0] - sorted_blocks_coords_in_reference[i][1])

            # B2.1 In the contig
            for i in range(0, len(sorted_blocks_coords_in_contig)-1):
                distance_between_blocks_contig.append(
                    sorted_blocks_coords_in_contig[i+1][0] - sorted_blocks_coords_in_contig[i][1])

            # C. Classification

            # C.1 Translocation - blocks over 1000 bps in distance in reference
            if any(i > 1000 for i in distances_between_blocks_ref):
                misassembly_list.append("translocation")

            # C.2 Insertion -
            if any(i > 50 for i in distance_between_blocks_contig):
                misassembly_list.append("insertion")

            # C.3 Deletion -
            if any(1000 >= i > 50 for i in distances_between_blocks_ref) and any(i <= 0 for i in distance_between_blocks_contig):
                misassembly_list.append("deletion")

            # C.4 - Duplication
            if any(i < 0 for i in distances_between_blocks_ref) and any(i <= 0 for i in distance_between_blocks_contig):
                misassembly_list.append("duplication")

            # C.5 - Rearangement
            if order_in_ref != sorted(order_in_ref):
                misassembly_list.append("rearrangement")

        # E - Catch all
        if len(misassembly_list) == 0:
            misassembly_list.append("inconsistency")

        missassembled_contigs[contig] = {'misassembly': misassembly_list,
                                         'contig length': contig_len,
                                         "n blocks": n_blocks,
                                         "distance_in_ref": distances_between_blocks_ref,
                                         "distance_in_contig": distance_between_blocks_contig,
                                         "reference": reference,
                                         "strands": strands,
                                         "blocks_coords_in_contig": blocks_coords_in_contig,
                                         "blocks_coords_in_reference": blocks_coords_in_reference,
                                         "order_in_ref": order_in_ref}

    return missassembled_contigs


def make_plot(mis_contigs, sample_id, assembler):
    """[summary]

    Parameters
    ----------
    mis_contigs : dictionary
        dictionary of misassembled contigs, with classification
    sample_id : string
        sample ID
    assembler : string
        assembler name
    """

    x = []  # contig lengths
    y = []  # n blocks
    z = []  # misassembly type
    w = []  # contig id
    for item, value in mis_contigs.items():
        x.append(value['contig length'])
        y.append(value['n blocks'])
        z.append(', '.join(value['misassembly']))
        w.append(item)

    df = pd.DataFrame(list(zip(x, y, z, w)),
                      columns=['Contig Length', 'n blocks', 'Misassembly', "Contig ID"])

    df['text'] = '<b>' + df['Misassembly'] + \
        '</b><br><br>Contig Name: ' + df['Contig ID'] + '<br>'

    trace = go.Scatter(x=df['Contig Length'],
                       y=df['n blocks'],
                       name=assembler, text=df['text'],
                       mode='markers',
                       opacity=0.7,
                       hovertemplate="%{text}" +
                       "Contig Length: %{x:.0f}bp<br>" +
                       "Fragments: %{y:.0}<br>" +
                       "<extra></extra>",)

    with open('{}_{}_trace.pkl'.format(sample_id, assembler), 'wb') as f:
        pickle.dump(trace, f)

    with open('{}_{}_contig_lenght.pkl'.format(sample_id, assembler), 'wb') as f:
        pickle.dump(x, f)


def make_df(sample_id, assembler, mis_contigs, filtered_paf_dict):
    """[summary]

    Parameters
    ----------
    sample_id : Str
        String with sample ID
    assembler : Str
        String with Assembler name
    mis_contigs : dictionary
        dictionary of misassembled contigs, with classification
    filtered_paf_dict: dictionary

    """

    df = pd.DataFrame(columns=["Contig", 'Sample', 'Reference',
                               'Ref Start', 'Ref End', 'Misassembly', "Assembler", 'Reference Length'])

    for contig_id in mis_contigs.keys():
        for contig_info in filtered_paf_dict[contig_id]:
            df = df.append({'Contig': contig_id, 'Sample': sample_id, 'Assembler': assembler, 'Reference': contig_info['reference'],
                            'Ref Start': contig_info['target start'], 'Ref End': contig_info['target end'],
                            'Misassembly': mis_contigs[contig_id]['misassembly'], 'Reference Length': contig_info['reference length']},
                           ignore_index=True)

    df.to_csv(sample_id + '_' + assembler + '_misassembly.csv')


def main(sample_id, assembler, assembly, mapping):

    # parse paf file
    paf_dict = parse_paf(mapping)

    # filter for contigs broken into multiple alignment blocks
    filtered_paf_dict = filter_dict(paf_dict)

    # classify misassembly
    mis_contigs = classify_misassembled_contigs(filtered_paf_dict)

    # global report
    mis_events = 0
    for contig in mis_contigs.keys():
        mis_events += len(mis_contigs[contig]['misassembly'])

    report_data = {"sample": sample_id, "assembler": assembler,
                   "misassembled_contigs": len(mis_contigs.keys()),
                   "misassembly_events": mis_events}
    
    logger.debug("Misassembly global report: {}".format(report_data))

    # reference report
    reference_report = {"sample": sample_id,
                        "assembler": assembler, "reference": {}}

    for contig in mis_contigs.keys():
        for reference in mis_contigs[contig]['reference']:
            if reference not in reference_report['reference'].keys():
                reference_report['reference'][reference] = {"misassembled_contigs": 1, "misassembly_events": len(mis_contigs[contig]['misassembly'])}
            else:
                reference_report['reference'][reference]["misassembled_contigs"] += 1
                reference_report['reference'][reference]["misassembly_events"] += len(mis_contigs[contig]['misassembly'])
    
    logger.debug("Misassembly reference report: {}".format(reference_report))

    # make plot trace - global
    make_plot(mis_contigs, sample_id, assembler)

    # prepare df for plot per reference
    make_df(sample_id, assembler, mis_contigs, filtered_paf_dict)

    # Write files for report
    with open("{}_{}_misassembly.json".format(sample_id, assembler), "w") as json_report:
        json_report.write(json.dumps(report_data, separators=(",", ":")))

    with open("{}_{}_misassembled_reference.json".format(sample_id, assembler), "w") as misassembly_dict:
        misassembly_dict.write(json.dumps(
            reference_report, separators=(",", ":")))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING)
