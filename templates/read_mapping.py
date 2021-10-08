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

import sys
import os
import subprocess
from subprocess import PIPE
import glob
import gzip
import json
import csv
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "03.11.2020"
__template__ = "READ_MAPPING-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    ASSEMBLER = '$assembler'
    ASSEMBLY = '$assembly'
    FILTERED_ASSEMBLY = '$filtered_assembly'
    FASTQ = '$params.fastq'
    BASEDIR = '$baseDir'
    THRESHOLD = '$params.mapped_reads_threshold'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("SAMPLE_ID: {}".format(SAMPLE_ID))
    logger.debug("ASSEMBLER: {}".format(ASSEMBLER))
    logger.debug("ASSEMBLY: {}".format(ASSEMBLY))
    logger.debug("FILTERED_ASSEMBLY: {}".format(FILTERED_ASSEMBLY))
    logger.debug("FASTQ: {}".format(FASTQ))
    logger.debug("BASEDIR: {}".format(BASEDIR))
    logger.debug("THRESHOLD: {}".format(THRESHOLD))

def map_to_assembly(assembly, reads, sample_id, assembler,threshold, n_reads_total):

    mapped_reads = 0

    cli = [
        "minimap2",
        "--sr",
        "-k21",
        "-N5"
        "--secondary=no",
        assembly,
        reads[0],
        reads[1],
        ">",
        "{}_{}_read_mapping_original.paf".format(sample_id, assembler)
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
    logger.info("Finished minimap2 with return code: {}".format(p.returncode))

    if p.returncode == 0:
        # count number of reads mapping
        with open("{}_{}_read_mapping_original.paf".format(sample_id, assembler)) as fh:
            csv_paf_file = csv.reader(fh, delimiter='\t')
            n_reads_mapping = 0
            for row in csv_paf_file:
                if int(row[10]) >= (int(row[1]) * float(threshold)):
                    n_reads_mapping += 1
        logger.debug("Number of reads mapping to assembly: {}".format(n_reads_mapping))


        try:
            mapped_reads = n_reads_mapping / n_reads_total
            logger.debug("Percentage of mapped reads to assembly: {}".format(mapped_reads))
        except ZeroDivisionError:
            mapped_reads = 0

    return mapped_reads


def main(sample_id, assembler, assembly, filtered_assembly, fastq, basedir, threshold):
    # get correct fastq files from directory
    all_readfiles = glob.glob(os.path.join(basedir, '/'.join(fastq.split('/')[:-1]), '*'))
    logger.debug("Read files found: {}".format(all_readfiles))
    reads = []
    for file in all_readfiles:
        if sample_id in file:
            reads.append(file)
    logger.debug("Matching read files: {}".format(reads))

    # get total number of reads
    n_reads_total = (sum(1 for line in gzip.open(reads[0], 'rb'))/4)+(sum(1 for line in gzip.open(reads[1], 'rb'))/4)
    logger.debug("Number of reads in fastq file: {}".format(n_reads_total))


    # map to original assembly
    mapped_reads_original = map_to_assembly(assembly, reads, sample_id, assembler,threshold, n_reads_total)
    with open("{}_{}_read_mapping_original.txt".format(sample_id, assembler), 'w') as fh:
        fh.write(str(mapped_reads_original * 100))

    # map to filtered assembly
    mapped_reads_filtered = map_to_assembly(filtered_assembly, reads, sample_id, assembler,threshold, n_reads_total)
    with open("{}_{}_read_mapping_filtered.txt".format(sample_id, assembler), 'w') as fh:
        fh.write(str(mapped_reads_filtered * 100))

    with open("{}_{}_read_mapping_report.json".format(sample_id, assembler), "w") as json_report:
        json_dic = {
            sample_id: {
                "assembler": assembler,
                "mapped_reads_original": mapped_reads_original * 100,
                "mapped_reads_filtered": mapped_reads_filtered * 100
        }}
        json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == '__main__':
    main(SAMPLE_ID, ASSEMBLER, ASSEMBLY, FILTERED_ASSEMBLY, FASTQ, BASEDIR, THRESHOLD)
