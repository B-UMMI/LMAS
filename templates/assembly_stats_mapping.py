#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to obtain information of percentage of mapping contigs/basepairs from a filtered assembly.
Collects information to produce a descriptive boxplot for mapped and unmapped contigs per assembler.

Metrics collected:
    * % mapped contigs - number of contigs that map to the reference files (and % over contigs > 1000bp)
    * % mapped bp - number of basepairs that map to the reference files (and % over total basepairs in
contigs > 1000bp)
    * Contiguity - largest % of reference covered by a single contig
    * lowest identity - % of identity to the reference of the worst mapping contig
    * breadth of coverage - % of the reference genome covered by the contigs
    * aligned contigs - number of aligned contigs to the reference
    * aligned basepairs - total number of basepairs aligned to te reference
    * C90 - Number of contigs, ordered by length, that cover 90% of the reference genome

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
from itertools import groupby
import utils

__version__ = "0.0.1"
__build__ = "28.10.2020"
__template__ = "ASSEMBLY_STATS_MAPPING-nf"

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


def get_covered_bases(covered_bases_list, ref_len):
    """
    Get ration of referee lengths (adjusted for triple reference) covered by mapping contigs
    :param covered_bases_list: list with alignment coordinates
    :param ref_len: expected reference length
    :return: % of reference covered by the alignment
    """
    sorted_list = sorted(covered_bases_list, key=lambda x: x[0])

    covered_bases = set()

    for item in sorted_list:
        start, stop = map(int, item[:])

        # Due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth.
        # Therefore, the coordinates are adjusted as follows:
        # [0; ref_len][ref_len+1; 2*ref_len][(2*ref_len)+1; 3*ref_len]
        for base in range(start, stop):
            if base <= ref_len:
                covered_bases.add(base)
            elif base <= 2*ref_len:
                covered_bases.add(base-ref_len)
            else:
                covered_bases.add(base-(2*ref_len))
    return len(covered_bases)/ref_len


def get_expanded_cigar(cigar):
    """

    :param cigar: string with cigar values
    :return: string with expanded cigar values for alignment
    """
    expanded_cigar = []
    cigar_parts = re.findall(r'\\d+[IDX=]', cigar)
    for cigar_part in cigar_parts:
        num = int(cigar_part[:-1])
        letter = cigar_part[-1]
        expanded_cigar.append(letter * num)
    return ''.join(expanded_cigar)


def get_lowest_window_identity(cigar, window_size):
    """

    :param cigar: string with alignment cigar
    :param window_size: int with window size
    :return: float with lowest identity value for mapping contigs
    """
    lowest_window_id = float('inf')
    expanded_cigar = get_expanded_cigar(cigar)

    for i in range(0, len(expanded_cigar) - window_size):
        cigar_window = expanded_cigar[i:i+window_size]
        window_id = cigar_window.count('=') / window_size
        if window_id < lowest_window_id:
            lowest_window_id = window_id
    if lowest_window_id == float('inf'):
        return 0.0
    return lowest_window_id


def get_phred_quality_score(identity):
    """
    Using the formula -log10(1-identity)*10, receives the identity of a contig and outputs the corresponding phred
    quality score.
    If the identity is 1, a Phred score of 60 (error rate of 0.0001%) is returned
    :param identity: float with identity values for the contig
    :return: float with phred score for contig identity
    """
    return - math.log10(1-identity) * 10 if identity < 1 else 60


def get_alignment_stats(paf_filename, ref_name, ref_length, df_phred):
    """
    Function to process the mapping (*.paf) file for a given reference.
    :param paf_filename: tabular file with alignment information for an assembler
    :param ref_name: reference name to filter from the paf_filename
    :param ref_length: expected reference length
    :param df_phred:
    :return:
        - contiguity: largest % of reference covered by a single contig
        - coverage:  % of the reference genome covered by the contigs (breadth of coverage)
        - lowest_identity: % of identity to the reference of the worst mapping contig
        - nID: Normalized identity by contig lenght
    """

    # Tracks the longest single alignment, in terms of the reference bases.
    covered_bases = []

    n_identity = []

    longest_alignment = 0

    alignment_dict = {'Reference': utils.REFERENCE_DIC[ref_name], 'Reference_Length': ref_length, 'Longest_Alignment': 0,
                     'Longest_Alignment_Cigar': '', 'Contigs': {}}

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                # parse values from PAF file
                contig_name, contig_length = parts[0], int(parts[1])
                start, end = int(parts[7]), int(parts[8])

                # number of residue matches, alignment block length
                matching_bases, total_bases = int(parts[9]), int(parts[10])
                cigar = parts[-1]

                if contig_name not in alignment_dict['Contigs'].keys():
                    alignment_dict['Contigs'][contig_name] = {'Length': contig_length, 'Base_Matches': matching_bases,
                                                             'Identity': None, 'Phred': None}
                else:
                    alignment_dict['Contigs'][contig_name]['Base_Matches'] += matching_bases

                if end - start > longest_alignment:
                    longest_alignment = end - start
                    longest_alignment_cigar = cigar
                longest_alignment = max(longest_alignment, end - start)
                covered_bases.append([start, end])

    # Calculate identity for all the contigs:
    for contig in alignment_dict['Contigs'].keys():
        alignment_dict['Contigs'][contig]['Identity'] = alignment_dict['Contigs'][contig]['Base_Matches'] / \
                                                       alignment_dict['Contigs'][contig]['Length']
        n_identity.append(alignment_dict['Contigs'][contig]['Base_Matches'])

        alignment_dict['Contigs'][contig]['Phred'] = get_phred_quality_score(alignment_dict['Contigs'][contig]['Identity'])
        df_phred = df_phred.append({'Assembler': os.path.basename(paf_filename).split('.')[0].rsplit('_')[-1],
                                    'Reference': alignment_dict['Reference'],
                                    'Contig': contig,
                                    'Contig Length': alignment_dict['Contigs'][contig]['Length'],
                                    'Phred Quality Score': alignment_dict['Contigs'][contig]['Phred']
                                    }, ignore_index=True)

    contiguity = longest_alignment / ref_length
    lowest_identity = get_lowest_window_identity(longest_alignment_cigar, 1000)

    coverage = get_covered_bases(covered_bases, ref_length)

    identity = sum(n_identity)/len(n_identity)

    return contiguity, coverage, lowest_identity, identity, df_phred


def get_c90(alignment_lengths, ref_len):
    """
    Returns the number of contigs, ordered by length, that cover at least 90% of the reference sequence.
    :param alignment_lengths: list with length of mapped contigs for the reference
    :param ref_len: int with the expected reference length
    :return: int with the number of contigs that represent
    """
    sorted_lengths = sorted(alignment_lengths, reverse=True)  # from longest to shortest
    target_length = ref_len * 0.9

    length_so_far = 0
    c90 = 0
    for contig_length in sorted_lengths:
        length_so_far += contig_length
        if length_so_far >= target_length:
            c90 += 1
    return c90


def parse_paf_files(sample_id, df, mapping, reference, assembler):
    """
    Parses fasta, paf files references and returns info in dataframe.
    :param sample_id: string with sample identifier
    :param df: pandas DataFrame with assembly stats
    :param mapping: paf file
    :param reference: path triple reference fasta file
    :param assembler: string with assembler name
    :return: pandas Dataframe with columns Reference, Assembler and C90
    """

    # Dataframe for C90 plot
    df_c90 = pd.DataFrame(columns=['Reference', 'Assembler', 'C90'])

    # Dataframe for Phred Score plot
    df_phred = pd.DataFrame(columns=['Assembler', 'Reference', 'Contig', 'Contig Length', 'Phred Quality Score'])

    # Mapping stats dict
    mapping_stats_dict = {'tableRow': [{
        "sample": sample_id,
        "assembler": assembler,
        "data": []
    }]}

    # filter dataframe for the assembler
    df_assembler = df[df['Assembler'] == assembler]

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(open(reference, "r"), lambda line: line[0] == ">"))

    fh = open(sample_id + '_' + assembler + "_breadth_of_coverage_contigs.csv", "w")
    fh.write("Reference, Breadth of Coverage, Contigs\\n")

    for header in references:
        header_str = header.__next__()[1:].strip().split()[0]
        reference_name = utils.REFERENCE_DIC[header_str]
        seq = "".join(s.strip() for s in references.__next__())

        df_assembler_reference = df_assembler[df_assembler['Mapped'] == header_str]

        mapped_contigs = df_assembler_reference['Contig Len'].astype('int').tolist()

        na50 = utils.get_N50(mapped_contigs)
        c90 = get_c90(mapped_contigs, len(seq)/3)  # adjust for triple reference
        df_c90 = df_c90.append({'Reference': reference_name, 'Assembler': assembler, 'C90': c90}, ignore_index=True)

        contiguity, coverage, lowest_identity, identity, df_phred = get_alignment_stats(mapping,
                                                                                        header_str,
                                                                                        len(seq)/3,
                                                                                        df_phred)

        fh.write(','.join([reference_name, str(coverage), str(len(mapped_contigs))]) + '\n')

        # Mapping stats dict
        mapping_ref_dict = {{'header': 'Reference',
                             'value': reference_name,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Reference length',
                             'value': len(seq)/3,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Contiguity',
                             'value': contiguity,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Identity',
                             'value': identity,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Lowest identity',
                             'value': lowest_identity,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Breadth of coverage',
                             'value': coverage,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'C90',
                             'value': c90,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Aligned contigs',
                             'value': len(mapped_contigs),
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'NA50',
                             'value': na50,
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id},
                            {'header': 'Aligned basepairs',
                             'value': sum(mapped_contigs),
                             'table': 'assembly_mapping_stats',
                             'sample': sample_id}
                            }

        mapping_stats_dict['data'].append(mapping_ref_dict)

        fh.close()

    return df_c90, df_phred, mapping_stats_dict


def main(sample_id, assembler, assembly, mapping, reference):

    # Dataframe with assembly info
    df = utils.parse_assemblies(assembler, assembly, mapping)

    to_plot_c90, to_plot_phred, json_dic = parse_paf_files(sample_id, df, mapping, reference, assembler)

    with open(".report.json", "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    to_plot_c90.to_csv(sample_id + '_' + assembler + '_c90.csv')
    to_plot_phred.to_csv(sample_id + '_' + assembler + '_phred.csv')

if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE)
