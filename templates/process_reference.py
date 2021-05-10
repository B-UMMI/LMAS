#!/usr/bin/env python3
"""
Purpose
-------
This script takes the multifasta file with the reference sequences and
converts them to tripled replicon reference sequences for assembly 
quality metric processing

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``reference_fasta``: path file to reference sequence
    - e.g.: ``'data/reference/*.fasta'``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
from itertools import groupby
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "07.05.2021"
__template__ = "PROCESS_REFERENCE-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    REFERENCE = '$reference_fasta'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("REFERENCE: {}".format(REFERENCE))

def main(reference):

    fasta_iter = utils.fasta_iter(reference)

    fh = open(reference)

    with open("triple_reference.fasta", "w") as fh:
        for header, seq in fasta_iter:
            print(header)
            fh.write('>' + header + '\\n')
            fh.write(seq * 3 + '\\n')

if __name__ == '__main__':
    main(REFERENCE)
