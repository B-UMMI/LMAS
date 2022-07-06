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
import json
import collections
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.2"
__build__ = "23.11.2020"
__template__ = "PROCESS_ASSEMBLY_STATS_GLOBAL-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    ASSEMBLY_STATS_GLOBAL_FILE = '$assembly_stats_global_files'.split()
    ASSEMBLY_STATS_GLOBAL_FILE_JSON = "$json_report".split()
    N_TARGET = float("$params.n_target")
    ASSEMBLER_SKIP = {"ABySS":json.loads("$params.abyss"), "GATBMiniaPipeline": json.loads("$params.gatb_minia"), 
                       "MetaHipMer2": json.loads("$params.metahipmer2"), "MINIA": json.loads("$params.minia"), "MEGAHIT": json.loads("$params.megahit"), 
                       "metaSPAdes": json.loads("$params.metaspades"), "Unicycler": json.loads("$params.unicycler"), "SPAdes": json.loads("$params.spades"),
                       "SKESA": json.loads("$params.skesa"), "VelvetOptimiser": json.loads("$params.velvetoptimiser"), "IDBA-UD": json.loads("$params.idba")}
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE: {}".format(
        ASSEMBLY_STATS_GLOBAL_FILE))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE_JSON: {}".format(
        ASSEMBLY_STATS_GLOBAL_FILE_JSON))
    logger.debug("N_TARGET: {}".format(
        N_TARGET))
    logger.debug("ASSEMBLER_STATUS: {}".format(
        ASSEMBLER_SKIP))


def main(assembly_stats_global_file, stats_json, n_target, assembler_skip):

    # Write JSON file report
    with open("global_assembly_stats.json", "w") as json_report:
        main_json = {}
        for data_report in stats_json:
            with open(data_report) as f:
                json_data = json.load(f)
                assembler = json_data["assembler"]
                sample_id = json_data["sample_id"]
                data_global = json_data["global"]
                data_filtered = json_data["filtered"]

            if sample_id not in main_json.keys():
                main_json[sample_id] = {"GlobalTable": [{
                    "assembler": assembler,
                    "original": data_global,
                    "filtered": data_filtered
                }]}
            else:
                main_json[sample_id]["GlobalTable"].append({
                    "assembler": assembler,
                    "original": data_global,
                    "filtered": data_filtered
                })

        # Add missing data for non-successful assemblies
        for assembler in utils.ASSEMBLER_NAMES:
            if assembler_skip[assembler]:
                for sample_id in main_json.keys():
                    if not any(d['assembler'] == assembler for d in main_json[sample_id]["GlobalTable"]):
                        main_json[sample_id]["GlobalTable"].append(
                            {
                                "assembler": assembler,
                                "original": {
                                    "contigs": 0,
                                    "basepairs": 0,
                                    "max_contig_size": 0,
                                    "N{}".format(int(n_target*100)): 0,
                                    "mapped_reads": 0,
                                    "Ns": 0},
                                "filtered": {
                                    "contigs": 0,
                                    "basepairs": 0,
                                    "N{}".format(int(n_target*100)): 0,
                                    "mapped_reads": 0,
                                    "Ns": 0}
                            })

        # sort final dictionary by assembler
        for sample in main_json.keys():
            main_json[sample]["GlobalTable"] = sorted(
                main_json[sample]["GlobalTable"], key=lambda i: i['assembler'].upper())

        json_report.write(json.dumps(main_json, separators=(",", ":")))


if __name__ == '__main__':
    main(ASSEMBLY_STATS_GLOBAL_FILE, ASSEMBLY_STATS_GLOBAL_FILE_JSON, N_TARGET, ASSEMBLER_SKIP)
