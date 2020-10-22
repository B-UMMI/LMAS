#!/usr/bin/env python3

"""
Purpose
-------
This module is intended parse the results of the process_assembly_stats_global.py for one or more
samples.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``assembly_stats_global_file`` : Path to output file.
    - e.g.: ``'assembly_stats_global.tsv'``

Generated output
----------------
None
"""

import os
import utils

__version__ = "0.0.1"
__build__ = "22.10.2020"
__template__ = "ASSEMBLY_STATS_GLOBAL-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    ASSEMBLY_STATS_GLOBAL_FILE = '$assembly_stats_global_files'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE: {}".format(ASSEMBLY_STATS_GLOBAL_FILE))


def main(assembly_stats_global_file):

    data = {}
    for file in assembly_stats_global_file:
        sample = file.split('_')[0]

        if sample not in data.keys():
            data[sample] = [{"data": file}]
        else:
            data[sample].append({"data": file})

    for sample in data.keys():
        with open(sample + 'csv') as csv_file:
            csv_file.write(','.join(['Assembler', 'Contigs', 'Basepairs', 'Max contig size', 'N50',
                                     'contigs>1000bp (%)', 'Basepairs in contigs>1000bp (%)', 'N50 in contigs>1000bp']))
            for item in data[sample]:
                csv_file.write(item['data'])


if __name__ == '__main__':
    main(ASSEMBLY_STATS_GLOBAL_FILE)
