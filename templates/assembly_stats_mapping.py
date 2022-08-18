#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information of percentage of mapping contigs/basepairs from a filtered assembly.
Collects information to produce a descriptive boxplot for mapped and unmapped contigs per assembler.

Metrics collected:
    * contiguity - largest % of reference covered by a single contig;
    * multiplicity - ratio of the length of alignable assembled sequence to covered sequence on the reference
    * validity - ratio of the length of the alignable assembled sequence to total basepairs in the aligned contigs;
    * parsimony - cost of the assembly (multiplicity over validity);
    * identity - ratio of identical basepairs in all aligned contigs to the reference;
    * lowest_identity -  % of identity to the reference of the worst mapping contig;
    * breadth_of_coverage -  % of the reference genome covered by the contigs;
    * L target - Number of contigs, ordered by length, that cover the target % of the reference genome;
    * aligned_contigs - number of aligned contigs to the reference;
    * NA target - length of aligned contigs of that length or longer that covers the target percentage of the total length of the reference;
    * NG target - length of aligned contigs of that length or longer that covers the target percentage of the sequence of the reference;
    * aligned_basepairs - total number of basepairs aligned to te reference;
    * Ns - Number of uncalled bases in the contigs that map to the reference. 

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
- ``reference``: fasta file of the complete triple reference genomes
    - e.g.: ``'reference.fasta' ``
- ``n_target``: target percetange, in float, for the NA and NG metrics
    - e.g.: ``'0.5' ``
- ``l_target``: target percetange, in float, for the L metric
    - e.g.: ``'0.9' ``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import math
import re
import json
import pandas as pd
import numpy as np
from itertools import groupby
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.2"
__build__ = "23.09.2021"
__template__ = "ASSEMBLY_STATS_MAPPING-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    MAPPING = '$mapping'
    REFERENCE = '$reference'
    N_TARGET = float("$params.n_target")
    L_TARGET = float("$params.l_target")
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("MAPPING: {}".format(MAPPING))
    logger.debug("REFERENCE: {}".format(REFERENCE))
    logger.debug("N_TARGET: {}".format(N_TARGET))
    logger.debug("L_TARGET: {}".format(L_TARGET))


def get_adjusted_basematches(contig_coords, base_matches):
    """
    Get the sum of the length of the overlaps, if they exist, when a contig is 
    broken into multiple alignment blocks. 
    :param contig_coords: list of lists containing start (0 based, closed) and stop (0 based, open) positions
    :param base_matches: current number of base matches as reported in the PAF file
    :returns: int with adjusted number of base matches, removing the overlaps
    """
    overlap_len = sum(utils.get_check_overlap(contig_coords))
    return base_matches - overlap_len


def get_covered_bases(covered_bases_list, ref_len):
    """
    Get ration of reference lengths (adjusted for triple reference) covered by mapping contigs
    :param covered_bases_list: list with alignment coordinates
    :param ref_len: expected reference length
    :return: % of reference covered by the alignment
    """
    sorted_list = sorted(covered_bases_list, key=lambda x: x[0])

    covered_bases = set()

    for item in sorted_list:
        start, stop = map(int, item[:])

        # Due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth.
        for base in range(start, stop):
            covered_bases.add(utils.adjust_reference_coord(base, ref_len))

    if len(covered_bases) == 0:
        return 0 
    else: 
        return len(covered_bases) - 1 #0 based index


def get_aligned_bases(alignment_coods):
    """
    Get number of bases that align to a reference in a contig.
    :param alignment_coods: list containing list with query start in the contig,
                            query end in the contig
    :return: length of bases in a contig that align to a reference
    """
    sorted_list = sorted(alignment_coods, key=lambda x: x[0])

    aligned_bases = set()
    for item in sorted_list:
        start, stop = map(int, item[:])
        for base in range(start, stop):
            aligned_bases.add(base)

    alignment_block_len = len(aligned_bases)
    return alignment_block_len


def get_identity(n_identity):
    """
    Calculates the average and lowest identity from a list of values (identity of each individual contig)
    :param n_identity: list with identity of each mapping contig
    :return float with average identity, float with lowest identity
    """
    identity = (sum(n_identity)/len(n_identity)
                ) if len(n_identity) > 0 else 0
    lowest_identity = min(n_identity) if len(n_identity) > 0 else 0
    return identity, lowest_identity


def get_phred_quality_score(identity):
    """
    Using the formula -log10(1-identity)*10, receives the identity of a contig and outputs the corresponding phred
    quality score.
    If the identity is 1, a Phred score of 60 (error rate of 0.0001%) is returned
    :param identity: float with identity values for the contig
    :return: float with phred score for contig identity
    """
    return - math.log10(1-identity) * 10 if identity < 1 else 60


def mapping_stats(sample_id, assembler, df, mapping_list, n_target, l_target):
    """
    Parses mapping and assembly data and returns mapping statistics
    :param sample_id: string with sample identifier
    :param df: pandas DataFrame with assembly stats
    :param mapping: paf file
    :param reference: path triple reference fasta file
    :param assembler: string with assembler name
    :param n_target: Target percentage reference length for NA and NG metrics
    :param l_target: Target percentage reference length for L metric
    :return: pandas Dataframe with columns Reference, Assembler and C90
    """

    # Dataframe for Phred Score plot
    df_phred = pd.DataFrame(columns=[
                            'Assembler', 'Reference', 'Contig', 'Contig Length', 'Phred Quality Score'])

    # Dataframes for assembly stats
    df_na = pd.DataFrame(
        columns=['Reference', 'Assembler', 'NAx', 'Basepairs'])
    df_ng = pd.DataFrame(
        columns=['Reference', 'Assembler', 'NGx', 'Basepairs'])
    df_lx = pd.DataFrame(columns=['Reference', 'Assembler', 'Lx', 'nContigs'])

    df_coverage = pd.DataFrame(
        columns=['Reference', 'Breadth of Coverage', 'Contigs'])

    # Mapping stats dict
    mapping_stats_dict = {
        "sample_id": sample_id,
        "ReferenceTables": {}}

    # filter dataframe for the assembler
    df_assembler = df[df['Assembler'] == assembler]

    # calculate stats per reference
    for alignment_dict in mapping_list:

        logger.debug("Calculating stats for {}".format(
            alignment_dict['Reference']))

        df_assembler_reference = df_assembler[df_assembler['Mapped']
                                              == alignment_dict['Reference']]

        mapped_contigs = df_assembler_reference['Contig Len'].unique().astype(
            'int').tolist()

        # Uncalled bases
        Ns = sum(df_assembler_reference['#N'].astype('int').tolist())

        # Contiguity
        for x in np.arange(0.0, 1.01, 0.01):

            x = round(x, 2) 

            # NAx
            nax = utils.get_Nx(mapped_contigs, x)
            df_na = df_na.append({'Reference': alignment_dict['Reference'], 'Assembler': assembler,
                                  'NAx': round(x*100), 'Basepairs': nax}, ignore_index=True)
            # NGx
            ngx = utils.get_NGx(
                mapped_contigs, alignment_dict['Reference_Length'], x)
            df_ng = df_ng.append({'Reference': alignment_dict['Reference'], 'Assembler': assembler,
                                  'NGx': round(x*100), 'Basepairs': ngx}, ignore_index=True)
            # Lx
            lx = utils.get_Lx(
                mapped_contigs, alignment_dict['Reference_Length'], x)
            df_lx = df_lx.append({'Reference': alignment_dict['Reference'], 'Assembler': assembler,
                                  'Lx': round(x*100), 'nContigs': lx}, ignore_index=True)

        na50 = utils.get_Nx(mapped_contigs, n_target)
        ng50 = utils.get_NGx(
            mapped_contigs, alignment_dict['Reference_Length'], n_target)
        l90 = utils.get_Lx(
            mapped_contigs, alignment_dict['Reference_Length'], l_target)

        # Calculate identity and Phred Score for all the contigs:
        sum_contig_length = 0
        alignment_block_list = []
        n_identity = []
        for contig in alignment_dict['Contigs'].keys():
            sum_contig_length += alignment_dict['Contigs'][contig]['Length']

            # Sometimes this value is > 1 when a contig is split into multiple alignment blocks
            # that overlap
            alignment_dict['Contigs'][contig]['Identity'] = alignment_dict['Contigs'][contig]['Base_Matches'] / \
                alignment_dict['Contigs'][contig]['Length']
            if alignment_dict['Contigs'][contig]['Identity'] <= 1: # sanity test
                n_identity.append(alignment_dict['Contigs'][contig]['Identity'])
            else: # in case of an overlap, consider identity as 1 as all bases match the reference
                n_identity.append(1)

            alignment_dict['Contigs'][contig]['Phred'] = get_phred_quality_score(
                alignment_dict['Contigs'][contig]['Identity'])
            df_phred = df_phred.append({'Assembler': assembler,
                                        'Reference': alignment_dict['Reference'],
                                        'Contig': contig,
                                        'Contig Length': alignment_dict['Contigs'][contig]['Length'],
                                        'Phred Quality Score': alignment_dict['Contigs'][contig]['Phred']
                                        }, ignore_index=True)

            # If a contig is broken into multiple blocks
            # the coords are adjusted when calculating the length
            if len(alignment_dict['Contigs'][contig]['Alignment_Blocks_Coords']) > 1:
                alignment_block_len = get_aligned_bases(
                    alignment_dict['Contigs'][contig]['Alignment_Blocks_Coords'])
                alignment_block_list.append(alignment_block_len)
            else:
                alignment_block_len = alignment_dict['Contigs'][contig]['Alignment_Blocks_Coords'][0][1] - \
                    alignment_dict['Contigs'][contig]['Alignment_Blocks_Coords'][0][0]
                alignment_block_list.append(alignment_block_len)

        identity, lowest_identity = get_identity(n_identity)

        # Contiguity
        contiguity = alignment_dict['Longest_Alignment'] / \
            alignment_dict['Reference_Length']

        # COMPASS Metrics
        sum_ci = get_covered_bases(
            alignment_dict['Covered_Bases'], alignment_dict['Reference_Length'])
        sum_ri = alignment_dict['Reference_Length']
        sum_ai = sum(alignment_block_list)
        sum_si = sum_contig_length

        coverage = sum_ci / sum_ri

        validity = sum_ai / sum_si if sum_si != 0 else 0

        multiplicity = sum_ai / sum_ci if sum_ci != 0 else 0

        parsimony = sum_si / sum_ci if sum_ci != 0 else 0

        # Update Coverage Dataframe
        df_coverage = df_coverage.append({'Reference': alignment_dict['Reference'],
                                          'Breadth of Coverage': coverage, 'Contigs': len(mapped_contigs)}, ignore_index=True)

        # Update Mapping stats dict
        mapping_stats_dict["ReferenceTables"][alignment_dict['Reference']] = {
            "assembler": assembler,
            "contiguity": contiguity,
            "multiplicity": multiplicity,
            "validity": validity,
            "parsimony": parsimony,
            "identity": identity,
            "lowest_identity": lowest_identity,
            "breadth_of_coverage": coverage,
            "L90": l90,
            "aligned_contigs": len(mapped_contigs),
            "NA50": na50,
            "NG50": ng50,
            "aligned_basepairs": sum_ci,
            "Ns": Ns,
        }

        logger.debug("  - Stats json: {}".format(
            mapping_stats_dict["ReferenceTables"][alignment_dict['Reference']]))

    return df_na, df_ng, df_lx, df_phred, df_coverage, mapping_stats_dict


def parse_paf_file(paf_filename, reference):
    """
    Function to process the mapping (*.paf) file and return it as 
    a dictionary containing the information of the aligning contigs
    per reference 

    :param paf_filename: path for tabular file with alignment information for an assembler, in PAF format
    :param reference: path triple reference fasta file

    :return:
        - alignment_dictitonary_list: list of dictionary containing the information 
                                      on the mapped contigs per reference
    """
    alignment_dictitonary_list = []

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(
        open(reference, "r"), lambda line: line[0] == ">"))

    logger.debug("Processing mapping information for...")

    for header in references:
        reference_name = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in references.__next__())
        reference_length = int(len(seq)/3)

        logger.debug(
            "  - {} with {} basepairs".format(reference_name, reference_length))

        alignment_dict = {'Reference': reference_name,
                          'Reference_Length': reference_length, 'Longest_Alignment': 0,
                          'Covered_Bases': [], 'Contigs': {}}

        with open(paf_filename) as paf:
            for line in paf:
                parts = line.strip().split('\t')
                if parts[5] == reference_name:
                    # parse values from PAF file
                    # contig info
                    contig_name, contig_length = parts[0], int(parts[1])
                    contig_start, contig_end = int(parts[2]), int(parts[3])

                    # coords in reference
                    ref_start, ref_end = int(parts[7]), int(parts[8])

                    # number of residue matches
                    matching_bases = int(parts[9])

                    # alignment block length - save coords in contig
                    alignment_block_list = [
                        int(parts[2]), int(parts[3])]

                    if contig_name not in alignment_dict['Contigs'].keys():
                        alignment_dict['Contigs'][contig_name] = {'Length': contig_length, 'Base_Matches': matching_bases,
                                                                  'Identity': None, 'Phred': None, 'Covered_Bases': [[ref_start, ref_end]],
                                                                  'Alignment_Blocks_Coords': [alignment_block_list],
                                                                  'Contig_Coords': [[contig_start, contig_end]]}
                    else:
                        alignment_dict['Contigs'][contig_name]['Contig_Coords'].append(
                            [contig_start, contig_end])
                        alignment_dict['Contigs'][contig_name]['Base_Matches'] += matching_bases
                        alignment_dict['Contigs'][contig_name]['Covered_Bases'].append([
                            ref_start, ref_end])

                    alignment_dict['Longest_Alignment'] = max(
                        alignment_dict['Longest_Alignment'], ref_end - ref_start)
                    alignment_dict['Covered_Bases'].append(
                        [ref_start, ref_end])
                    alignment_dict['Contigs'][contig_name]['Alignment_Blocks_Coords'].append(
                        alignment_block_list)

        # calculate base matches compensating for possible overlap
        for contig_name in alignment_dict['Contigs']:
            if utils.check_overlap(alignment_dict['Contigs'][contig_name]['Contig_Coords']):
                adjusted_matching_bases = get_adjusted_basematches(
                    alignment_dict['Contigs'][contig_name]['Contig_Coords'], alignment_dict['Contigs'][contig_name]['Base_Matches'])
                alignment_dict['Contigs'][contig_name]['Base_Matches'] = adjusted_matching_bases

        alignment_dictitonary_list.append(alignment_dict)

    return alignment_dictitonary_list


def main(sample_id, assembler, assembly, mapping, reference, n_target, l_target):

    # Dataframe with assembly info
    df = utils.parse_assemblies(sample_id, assembler, assembly, mapping)
    df.to_csv(sample_id + '_' + assembler + '_df.csv')

    # list of dictionaries with mapping info per reference
    alignment_dictitonary_list = parse_paf_file(mapping, reference)

    # get mapping stats
    to_plot_nax, to_plot_ngx, to_plot_lx, to_plot_phred, to_plot_coverage, json_dic = mapping_stats(
        sample_id, assembler, df, alignment_dictitonary_list, n_target, l_target)

    # save output files
    with open("{}_{}_report.json".format(sample_id, assembler), "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    to_plot_nax.to_csv(sample_id + '_' + assembler + '_nax.csv')
    to_plot_ngx.to_csv(sample_id + '_' + assembler + '_ngx.csv')
    to_plot_lx.to_csv(sample_id + '_' + assembler + '_lx.csv')
    to_plot_phred.to_csv(sample_id + '_' + assembler + '_phred.csv')
    to_plot_coverage.to_csv(sample_id + '_' + assembler +
                            "_breadth_of_coverage_contigs.csv")


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE, N_TARGET, L_TARGET)
