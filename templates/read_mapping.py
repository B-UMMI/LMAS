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
from subprocess import PIPE
import glob
import gzip
import json
import utils

__version__ = "0.0.1"
__build__ = "03.11.2020"
__template__ = "READ_MAPPING-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    FASTQ = '$params.fastq'
    BASEDIR = '$baseDir'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("FASTQ: {}".format(FASTQ))
    logger.debug("BASEDIR: {}".format(BASEDIR))


def main(sample_id, assembler, assembly, fastq, basedir):
    # get correct fastq files from directory
    all_readfiles = glob.glob(os.path.join(basedir, '/'.join(fastq.split('/')[:-1]), '*'))
    logger.debug("Read files found: {}".format(all_readfiles))
    reads = []
    for file in all_readfiles:
        if sample_id in file:
            reads.append(file)
    logger.debug("Matching read files: {}".format(reads))

    # call minimap2
    cli = [
        "minimap2",
        "-x",
        "sr",
        assembly,
        reads[0],
        reads[1],
        '>',
        "{}_{}_read_mapping.paf".format(sample_id, assembler)
    ]

    logger.debug("Running minimap2 subprocess with command: {}".format(' '.join(cli)))

    p = subprocess.Popen(' '.join(cli), shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()

    try:
        stderr = stderr.decode("utf8")
        stdout = stdout.decode("utf8")
    except (UnicodeDecodeError, AttributeError):
        stderr = str(stderr)
        stdout = str(stdout)

    logger.info("Finished minimap2 subprocess with STDOUT:\\n"
                "======================================\\n{}".format(stdout))
    logger.info("Fished minimap2 subprocesswith STDERR:\\n"
                "======================================\\n{}".format(stderr))
    logger.info("Finished minimap2 with return code: {}".format(
        p.returncode))

    if p.returncode == 0:
        # count number of reads mapping
        n_reads_mapping = sum(1 for line in open("{}_{}_read_mapping.paf".format(sample_id, assembler)))

        # get number of reads
        n_reads = sum(1 for line in gzip.open(reads[0], 'rb'))/4

        with open("{}_{}_read_mapping.txt".format(sample_id, assembler), 'w') as fh:
            try:
                mapped_reads = n_reads_mapping/(n_reads * 2)
            except ZeroDivisionError:
                mapped_reads = 0
            fh.write(str(mapped_reads))

        with open(".report.json", "w") as json_report:
            json_dic = {
                "tableRow": [{
                    "sample": sample_id,
                    "data": [
                        {"header": "Mapped reads (%)",
                         "value": mapped_reads * 100,
                         "table": "assembly_global_stats",
                         "assembler": assembler},
                    ]
                }]
            }
            json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, FASTQ, BASEDIR)
