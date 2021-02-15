#!/usr/bin/env python3
"""
Purpose
-------
This script count reads in fastq files and return json with info for report
Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``sample_id``: Sample Identification string.
    - e.g.: ``'SampleA'``
- ``fastq``: glob path file to read files for matching with sample_id
    - e.g.: ``'data/fastq/*_{1,2}.*'``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import gzip
import json
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "15.02.2021"
__template__ = "PROCESS_READS-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    FASTQ = '$params.fastq'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("FASTQ: {}".format(FASTQ))


def main(sample_id, fastq):

    # get total number of reads
    n_reads_total = (sum(1 for line in gzip.open(fastq[0], 'rb'))/4)+(sum(1 for line in gzip.open(fastq[1], 'rb'))/4)

    with open("{}_reads_report.json".format(sample_id), "w") as json_report:
        json_dic = {
            sample_id: {
                "reads": n_reads_total}}
        json_report.write(json.dumps(json_dic, separators=(",", ":")))

if __name__ == '__main__':
    main(SAMPLE_ID, FASTQ)
