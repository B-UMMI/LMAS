#!/usr/bin/env python3
"""
Purpose
-------
This script reads the version of the assembler and saves them in an json
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``versions``: List with file paths with assembler versions as string
    - e.g.: ``[.a_version, .b_versions]``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import json
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "07.03.2021"
__template__ = "PROCESS_VERSION-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    VERSIONS = '$version'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("VERSIONS: {}".format(VERSIONS))


ASSEMBLER_PROCESS_DICT = {
    "ABySS": "ABYSS",
    "BCALM2": "BCALM2",
    "GATBMiniaPipeline": "GATBMINIAPIPELINE",
    "MetaHipMer2": "METAHIPMER2",
    "MINIA": "MINIA",
    "MEGAHIT": "MEGAHIT",
    'metaSPAdes': "METASPADES",
    "Unicycler": "UNICYCLER",
    "SPAdes": "SPADES",
    "SKESA": "SKESA",
    "VelvetOptimiser": "VELVETOPTIMISER",
    "IDBA": "IDBA"}


def main(versions):

    version_dict = {}
    for v_file in versions:
        assembler_name = ASSEMBLER_PROCESS_DICT[v_file.replace(
            '.', '').split('_')[-2]]
        with open(v_file) as fh:
            assembler_version = fh.readline().strip()
            if assembler_name not in version_dict.keys():
                version_dict[assembler_name] = assembler_version

    with open("versions.json", "w") as json_report:
        json_report.write(json.dumps(version_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(VERSIONS)
