#!/usr/bin/env python3
"""
Purpose
-------
This script is a wrapper around minimap2 to map the readfiles to the obtained assembly using minimap2

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
- ``fastq``: glob path file to read files for matching with sample_id
    - e.g.: ``'data/fastq/*_{1,2}.*'``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import subprocess
import utils

__version__ = "0.0.1"
__build__ = "03.11.2020"
__template__ = "READ_MAPPING-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    REFERENCE = '$reference'
    FASTQ = '$params.fastq'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("REFERENCE: {}".format(REFERENCE))
    logger.debug("FASTQ: {}".format(FASTQ))


def main(sample_id, assembler, assembly, reference, fastq):
    print(fastq)


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, REFERENCE, FASTQ)