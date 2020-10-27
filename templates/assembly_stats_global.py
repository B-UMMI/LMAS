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
import utils

__version__ = "0.0.1"
__build__ = "22.10.2020"
__template__ = "ASSEMBLY_STATS_GLOBAL-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))


def get_contig_lists(fasta):
    """
    From a fasta iterator, get lists with contig lengths
    :param fasta: yield tuples of header, sequence
    :return:
        - contig_len: list with all contig lenghts in the assembly
        - contigs_len_over_1000: list with contig lenghts filtered for a minimum of 1000 nucleotides
    """
    contigs_len_over_1000 = []
    contigs_len = []

    for header, seq in fasta:
        if len(seq) > 1000:
            contigs_len_over_1000.append(len(seq))
        contigs_len.append(len(seq))

    return contigs_len, contigs_len_over_1000


def main(sample_id, assembler, assembly):

    contigs, contigs_over_1000bp = get_contig_lists(utils.fasta_iter(assembly))

    n50_contigs = utils.get_N50(contigs)
    n50_contigs_over_1000bp = utils.get_N50(contigs_over_1000bp)

    with open(".report.json", "w") as json_report:
        json_dic = {
            "tableRow": [{
                "assembler": assembler,
                "data": [
                    {"header": "Contigs",
                     "value": len(contigs),
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "Basepairs",
                     "value": sum(contigs),
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "Max contig size",
                     "value": max(contigs) if len(contigs) > 0 else 0,
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "N50",
                     "value": n50_contigs,
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "contigs>1000bp (%)",
                     "value": len(contigs_over_1000bp),
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "Basepairs in contigs>1000bp (%)",
                     "value": sum(contigs_over_1000bp),
                     "table": "assembly_global_stats",
                     "sample": sample_id},
                    {"header": "N50 in contigs>1000bp",
                     "value": n50_contigs_over_1000bp,
                     "table": "assembly_global_stats",
                     "sample": sample_id}
                ]
            }]
        }

        json_report.write(json.dumps(json_dic, separators=(",", ":")))

    with open(sample_id + '_' + assembler + "_global_assembly_stats_global.csv", "w") as cvs_file:
        cvs_file.write(','.join([assembler, f'{len(contigs)}', f'{sum(contigs)}',
                                 f'{max(contigs)if len(contigs) > 0 else 0 }', f'{n50_contigs}',
                                 f'{len(contigs_over_1000bp)}', f'{sum(contigs_over_1000bp)}',
                                 f'{n50_contigs_over_1000bp}']))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY)
