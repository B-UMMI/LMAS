#!/usr/bin/env python3

"""
Purpose
-------
This module is intended parse the results of the process_assembly_stats_mapping.py for one or more
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
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "24.11.2020"
__template__ = "PROCESS_ASSEMBLY_STATS_MAPPINGL-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    ASSEMBLY_STATS_GLOBAL_FILE_JSON = "$json_report".split()
    ASSEMBLER_SKIP = {"ABySS":json.loads("$params.abyss"), "BCALM2": json.loads("$params.bcalm"), "GATBMiniaPipeline": json.loads("$params.gatb_minia"), 
                       "MetaHipMer2": json.loads("$params.metahipmer2"), "MINIA": json.loads("$params.minia"), "MEGAHIT": json.loads("$params.megahit"), 
                       "metaSPAdes": json.loads("$params.metaspades"), "Unicycler": json.loads("$params.unicycler"), "SPAdes": json.loads("$params.spades"),
                       "SKESA": json.loads("$params.skesa"), "VelvetOptimiser": json.loads("$params.velvetoptimiser"), "IDBA-UD": json.loads("$params.idba"),
                       "RAVEN": json.loads("$params.raven"),"FLYE": json.loads("$params.flye"), "METAFLYE": json.loads("$params.metaflye"),
                       "RA": json.loads("$params.ra"),"WTDBG2": json.loads("$params.wtdbg2")}
    MODE = "$params.wf"

    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE_JSON: {}".format(
        ASSEMBLY_STATS_GLOBAL_FILE_JSON))
    logger.debug("ASSEMBLER_SKIP: {}".format(
        ASSEMBLER_SKIP))
    logger.debug("MODE: {}".format(
        MODE))
    


def main(stats_json, assembler_skip, mode):

    if mode == "default" or mode == "Illumina" or mode == "illumina":
        assembler_skip["RAVEN"] = False
        assembler_skip["FLYE"] = False
        assembler_skip["METAFLYE"] = False
        assembler_skip["RA"] = False
        assembler_skip["WTDBG2"] = False

    elif mode == "ONT" or mode == "ont":
        assembler_skip["ABySS"] = False
        assembler_skip["BCALM2"] = False
        assembler_skip["GATBMiniaPipeline"] = False
        assembler_skip["MetaHipMer2"] = False
        assembler_skip["MINIA"] = False
        assembler_skip["MEGAHIT"] = False
        assembler_skip["metaSPAdes"] = False
        assembler_skip["Unicycler"] = False
        assembler_skip["SPAdes"] = False
        assembler_skip["SKESA"] = False
        assembler_skip["VelvetOptimiser"] = False
        assembler_skip["IDBA-UD"] = False
    
    elif mode == "Hybrid" or mode == "hybrid":
        assembler_skip["RAVEN"] = False
        assembler_skip["FLYE"] = False
        assembler_skip["METAFLYE"] = False
        assembler_skip["RA"] = False
        assembler_skip["WTDBG2"] = False
        assembler_skip["ABySS"] = False
        assembler_skip["BCALM2"] = False
        assembler_skip["GATBMiniaPipeline"] = False
        assembler_skip["MetaHipMer2"] = False
        assembler_skip["MINIA"] = False
        assembler_skip["MEGAHIT"] = False
        assembler_skip["SKESA"] = False
        assembler_skip["VelvetOptimiser"] = False
        assembler_skip["IDBA-UD"] = False

    else:
        pass

    # Write JSON file report
    with open("global_assembly_mapping_stats.json", "w") as json_report:
        main_json = {}

        for data_report in stats_json:
            with open(data_report) as f:
                json_data = json.load(f)
                sample_id = json_data["sample_id"]
                reference_table_data = json_data["ReferenceTables"]

            for reference_id, reference_mapping_stats in reference_table_data.items():
                if sample_id not in main_json.keys():
                    main_json[sample_id] = {"ReferenceTable": {
                        reference_id: [reference_mapping_stats]}}
                else:
                    if reference_id not in main_json[sample_id]["ReferenceTable"].keys():
                        main_json[sample_id]["ReferenceTable"][reference_id] = [
                            reference_mapping_stats]
                    else:
                        main_json[sample_id]["ReferenceTable"][reference_id].append(
                            reference_mapping_stats)

        # Add missing data for non-successful assemblies
        for assembler in utils.ASSEMBLER_NAMES:
            if assembler_skip[assembler]:
                for sample_id in main_json.keys():
                    for reference in main_json[sample_id]["ReferenceTable"]:
                        if not any(d['assembler'] == assembler for d in main_json[sample_id]["ReferenceTable"][reference]):
                            main_json[sample_id]["ReferenceTable"][reference].append(
                                {
                                    "assembler": assembler,
                                    "contiguity": 0,
                                    "multiplicity": 0,
                                    "validity": 0,
                                    "parsimony": 0,
                                    "identity": 0,
                                    "lowest_identity": 0,
                                    "breadth_of_coverage": 0,
                                    "L90": 0,
                                    "aligned_contigs": 0,
                                    "NA50": 0,
                                    "NG50": 0,
                                    "aligned_basepairs": 0,
                                    "Ns": 0,
                                })
        # sort final dictionary by assembler
        for sample in main_json.keys():
            for reference in main_json[sample]["ReferenceTable"]:
                main_json[sample]["ReferenceTable"][reference] = sorted(
                    main_json[sample]["ReferenceTable"][reference], key=lambda i: i['assembler'].upper())

        json_report.write(json.dumps(main_json, separators=(",", ":")))


if __name__ == '__main__':
    main(ASSEMBLY_STATS_GLOBAL_FILE_JSON, ASSEMBLER_SKIP, MODE)
