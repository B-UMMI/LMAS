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
import numpy as np
from itertools import groupby
try:
    import utils
except ImportError:
    from templates import utils

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

    return len(covered_bases)/ref_len, len(covered_bases)


def get_phred_quality_score(identity):
    """
    Using the formula -log10(1-identity)*10, receives the identity of a contig and outputs the corresponding phred
    quality score.
    If the identity is 1, a Phred score of 60 (error rate of 0.0001%) is returned
    :param identity: float with identity values for the contig
    :return: float with phred score for contig identity
    """
    return - math.log10(1-identity) * 10 if identity < 1 else 60


def get_multiplicity(covered_bases_list, ref_len):
    """
    Get ratio of sum of aligned blocks with total length of aligned blocks
    :param covered_bases_list: list with alignment coordinates
    :param ref_len: expected reference length
    :return: multiplicity
    """
    sorted_list = sorted(covered_bases_list, key=lambda x: x[0])

    covered_bases = set()
    total_bases = 0

    for item in sorted_list:
        start, stop = map(int, item[:])

        # Due to the triple reference, the values need to be adjusted as not to over-estimate coverage breadth.
        for base in range(start, stop):
            covered_bases.add(utils.adjust_reference_coord(base, ref_len))
            total_bases += 1
    if total_bases > 0:
        return total_bases / len(covered_bases)
    else:
        return 0


def get_validity(covered_bases_list, sum_contig_length):
    """
    Get ratio of sum of bases that map to the reference to sum of lenght of aligned contigs 
    :param covered_bases_list: list with alignment coordinates
    :param sum_contig_length: int of sum of total lenght of all aligned contigs
    :return: validity
    """
    sorted_list = sorted(covered_bases_list, key=lambda x: x[0])
    total_bases = 0

    for item in sorted_list:
        start, stop = map(int, item[:])
        for base in range(start, stop):
            total_bases += 1

    if total_bases > 0:
        return  total_bases / sum_contig_length
    else:
        return 0


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
    alignment_dict = {'Reference': utils.REFERENCE_DIC[ref_name], 'Reference_Length': ref_length,
                      'Longest_Alignment': 0, 'Contigs': {}}

    with open(paf_filename) as paf:
        for line in paf:
            parts = line.strip().split('\t')
            if parts[5] == ref_name:
                # parse values from PAF file
                contig_name, contig_length = parts[0], int(parts[1])
                start, end = int(parts[7]), int(parts[8])

                # number of residue matches
                matching_bases = int(parts[9])

                if contig_name not in alignment_dict['Contigs'].keys():
                    alignment_dict['Contigs'][contig_name] = {'Length': contig_length, 'Base_Matches': matching_bases,
                                                              'Identity': None, 'Phred': None}
                else:
                    alignment_dict['Contigs'][contig_name]['Base_Matches'] += matching_bases

                if end - start > longest_alignment:
                    longest_alignment = end - start
                longest_alignment = max(longest_alignment, end - start)
                covered_bases.append([start, end])

    # Calculate identity for all the contigs:
    sum_contig_length = 0
    for contig in alignment_dict['Contigs'].keys():
        sum_contig_length += alignment_dict['Contigs'][contig]['Length']

        alignment_dict['Contigs'][contig]['Identity'] = alignment_dict['Contigs'][contig]['Base_Matches'] / \
            alignment_dict['Contigs'][contig]['Length']
        n_identity.append(alignment_dict['Contigs'][contig]['Identity'])

        alignment_dict['Contigs'][contig]['Phred'] = get_phred_quality_score(
            alignment_dict['Contigs'][contig]['Identity'])
        df_phred = df_phred.append({'Assembler': os.path.basename(paf_filename).split('.')[0].rsplit('_')[-1],
                                    'Reference': alignment_dict['Reference'],
                                    'Contig': contig,
                                    'Contig Length': alignment_dict['Contigs'][contig]['Length'],
                                    'Phred Quality Score': alignment_dict['Contigs'][contig]['Phred']
                                    }, ignore_index=True)

    contiguity = longest_alignment / ref_length

    # COMPASS Metrics
    coverage, len_convered_bases = get_covered_bases(covered_bases, ref_length)
    print("coverage: {}".format(coverage))

    multiplicity = get_multiplicity(covered_bases, ref_length)
    print("multiplicity: {}".format(multiplicity))

    validity = get_validity(covered_bases, sum_contig_length)
    print("validity: {}".format(validity))

    parsimony = multiplicity / validity if validity != 0 else 0
    print("parsimony: {}".format(parsimony))

    identity = (sum(n_identity)/len(n_identity)) if len(n_identity) > 0 else 0
    lowest_identity = min(n_identity) if len(n_identity) > 0 else 0

    return contiguity, coverage, multiplicity, validity, parsimony, lowest_identity, identity, df_phred, len_convered_bases


def parse_paf_files(sample_id, df, mapping, reference, assembler, n_target, l_target):
    """
    Parses fasta, paf files references and returns info in dataframe.
    :param sample_id: string with sample identifier
    :param df: pandas DataFrame with assembly stats
    :param mapping: paf file
    :param reference: path triple reference fasta file
    :param assembler: string with assembler name
    :param n_target: Target percentage reference lenght for NA and NG metrics
    :param l_target: Target percentage reference lenght for L metric
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

    # Mapping stats dict
    mapping_stats_dict = {
        "sample_id": sample_id,
        "ReferenceTables": {}}

    # filter dataframe for the assembler
    df_assembler = df[df['Assembler'] == assembler]

    # iterator for reference files (sequence length is needed)
    references = (x[1] for x in groupby(
        open(reference, "r"), lambda line: line[0] == ">"))

    fh = open(sample_id + '_' + assembler +
              "_breadth_of_coverage_contigs.csv", "w")
    fh.write("Reference,Breadth of Coverage,Contigs\\n")

    for header in references:
        header_str = header.__next__()[1:].strip().split()[0]
        reference_name = utils.REFERENCE_DIC[header_str]
        seq = "".join(s.strip() for s in references.__next__())

        df_assembler_reference = df_assembler[df_assembler['Mapped'] == header_str]

        mapped_contigs = df_assembler_reference['Contig Len'].astype(
            'int').tolist()

        Ns = sum(df_assembler_reference['#N'].astype('int').tolist())

        # Assembly metrics
        for x in np.linspace(0, 1, 10):
            # NAx
            nax = utils.get_Nx(mapped_contigs, x)
            df_na = df_na.append({'Reference': reference_name, 'Assembler': assembler,
                                  'NAx': x, 'Basepairs': nax}, ignore_index=True)
            # NGx
            ngx = utils.get_NGx(mapped_contigs, len(seq)/3, x)
            df_ng = df_ng.append({'Reference': reference_name, 'Assembler': assembler,
                                  'NGx': x, 'Basepairs': ngx}, ignore_index=True)
            # Lx
            lx = utils.get_Lx(mapped_contigs, len(seq)/3, x)
            df_lx = df_lx.append({'Reference': reference_name, 'Assembler': assembler,
                                  'Lx': x, 'nContigs': lx}, ignore_index=True)

        na50 = utils.get_Nx(mapped_contigs, n_target)
        ng50 = utils.get_NGx(mapped_contigs, len(seq)/3, n_target)
        l90 = utils.get_Lx(mapped_contigs, len(seq)/3, l_target)

        contiguity, coverage, multiplicity, validity, parsimony, lowest_identity, identity, df_phred, covered_bases = get_alignment_stats(mapping,
                                                                                                                                          header_str,
                                                                                                                                          len(
                                                                                                                                              seq)/3,
                                                                                                                                          df_phred)

        fh.write(','.join([reference_name, str(coverage),
                           str(len(mapped_contigs))]) + '\\n')

        # Mapping stats dict
        mapping_stats_dict["ReferenceTables"][reference_name] = {
            "assembler": assembler,
            "contiguity": contiguity,
            "multiplicity": multiplicity,
            "validity": validity,
            "parsimony": parsimony,
            "identity": identity,
            "lowest_identity": lowest_identity,
            "breadth_of_coverage": coverage,
            "L90": l90 if l90 is not None else 0,
            "aligned_contigs": len(mapped_contigs),
            "NA50": na50,
            "NG50": ng50,
            "aligned_basepairs": covered_bases,
            "Ns": Ns,
        }

    fh.close()

    return df_na, df_ng, df_lx, df_phred, mapping_stats_dict


def main(sample_id, assembler, assembly, mapping, reference, n_target, l_target):

    # Dataframe with assembly info
    df = utils.parse_assemblies(sample_id, assembler, assembly, mapping)

    # save dataframe
    df.to_csv(sample_id + '_' + assembler + '_df.csv')

    to_plot_nax, to_plot_ngx, to_plot_lx, to_plot_phred, json_dic = parse_paf_files(sample_id, df,
                                                                                    mapping, reference, assembler,
                                                                                    n_target, l_target)

    with open("{}_{}_report.json".format(sample_id, assembler), "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    to_plot_nax.to_csv(sample_id + '_' + assembler + '_nax.csv')
    to_plot_ngx.to_csv(sample_id + '_' + assembler + '_ngx.csv')
    to_plot_lx.to_csv(sample_id + '_' + assembler + '_lx.csv')
    to_plot_phred.to_csv(sample_id + '_' + assembler + '_phred.csv')


if __name__ == '__main__':
    #main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, MAPPING, REFERENCE, N_TARGET, L_TARGET)
    main("mockSample", "GATBMiniaPipeline", "filtered_ERR2935805_GATBMiniaPipeline.fasta", "ERR2935805_GATBMiniaPipeline.paf",
    "Zymos_Genomes_triple_chromosomes.fasta", 0.5, 0.9)
