#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Purpose
-------
Get the basic statistics from a RAW assembly.
Outputs the statistics for all the contigs in the assembly and for the contigs with over 1000pb into a tsv and json file.

Assembly metrics collected: 
  * Assembler - assembler name
  * Contigs - number of contigs in the assembly
  * basepairs - total number of nucleotides in the assembly
  * Max contig size - size of the largest contig in the assembly
  * n50 - sequence length of the shortest contig at 50% of the total assembly length
  * contigs>1000bp (%) - number of contigs with size >= 1000 bp (and % over "Contigs")
  * bp in contigs>1000bp (%) - total number of nucleotides in contigs with size >= 1000 bp (and % over "basepairs")
  * n50 in contigs>1000bp - sequence length of the shortest contig at 50% of the total length of contigs with ´
size >= 1000 bp

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``assembler``: String with assembler name.
    - e.g.: ``'SPAdes'``
- ``assembly``: fasta file from the assembler
    - e.g.: ``'spades.fasta'``
- ``min_len``: value for minimum contig length
    - e.g.: ``'1000'``

Generated output
----------------
- ``.report.json``: Data structure for the report

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import json
import re
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "22.10.2020"
__template__ = "ASSEMBLY_STATS_GLOBAL-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    READ_MAPPING_STATS = '$read_mapping'
    MIN_LEN = "$params.minLength"
    N_TARGET = float("$params.n_target")
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("READ_MAPPING_STATS: {}".format(READ_MAPPING_STATS))
    logger.debug("MIN_LEN: {}".format(MIN_LEN))
    logger.debug("N_TARGET: {}".format(N_TARGET))


def get_contig_lists(fasta, min_len):
    """
    From a fasta iterator, get lists with contig lengths
    :param fasta: yield tuples of header, sequence
    :param min_len: minimum contig lenght
    :return:
        - contig_len: list with all contig lenghts in the assembly
        - contigs_len_over_1000: list with contig lenghts filtered for a minimum of 1000 nucleotides
    """
    contigs_len_over_1000 = []
    contigs_len = []
    Ns_all = 0
    Ns_over_1000 = 0

    for header, seq in fasta:
        Ns_all += len(re.findall("N", seq.upper()))
        if len(seq) > int(min_len):
            contigs_len_over_1000.append(len(seq))
            Ns_over_1000 += len(re.findall("N", seq.upper()))
        contigs_len.append(len(seq))

    return contigs_len, contigs_len_over_1000, Ns_all, Ns_over_1000


def main(sample_id, assembler, assembly, read_mapping_stats, min_len, n_target):
    

    contigs, contigs_over_min_len, Ns_all, Ns_over_1000 = get_contig_lists(utils.fasta_iter(assembly), min_len)

    n50_contigs = utils.get_Nx(contigs, n_target)
    n50_contigs_over_min_len = utils.get_Nx(contigs_over_min_len, n_target)

    # get read mapping stats
    with open(read_mapping_stats) as f:
        assembly_stats_json = json.load(f)
        if assembly_stats_json[sample_id]["assembler"] == assembler:
            mapped_reads_all = assembly_stats_json[sample_id]["mapped_reads_original"]
            mapped_reads_filtered = assembly_stats_json[sample_id]["mapped_reads_filtered"]
        else:
            mapped_reads_all = 0
            mapped_reads_filtered = 0
            logger.error(assembly_stats_json)

    with open("{}_{}_report.json".format(sample_id, assembler), "w") as json_report:
        json_dic = {
                "assembler": assembler,
                "sample_id": sample_id,
                "global": {
                    "contigs": len(contigs),
                    "basepairs": sum(contigs),
                    "max_contig_size": max(contigs) if len(contigs) > 0 else 0,
                    "N{}".format(int(n_target*100)): n50_contigs,
                    "mapped_reads": mapped_reads_all,
                    "Ns": Ns_all},
                "filtered": {
                        "contigs": len(contigs_over_min_len),
                        "basepairs": sum(contigs_over_min_len),
                        "N{}".format(int(n_target*100)): n50_contigs_over_min_len,
                        "mapped_reads": mapped_reads_filtered,
                        "Ns": Ns_over_1000}
                }

        json_report.write(json.dumps(json_dic, separators=(",", ":")))

        with open(sample_id + '_' + assembler + "_global_assembly_stats_global.csv", "w") as cvs_file:
            cvs_file.write(','.join([assembler, f'{len(contigs)}', f'{sum(contigs)}',
                                    f'{max(contigs)if len(contigs) > 0 else 0 }', f'{n50_contigs}',
                                    f'{len(contigs_over_min_len)}', f'{sum(contigs_over_min_len)}',
                                    f'{n50_contigs_over_min_len}']))

if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY,READ_MAPPING_STATS, MIN_LEN, N_TARGET)
