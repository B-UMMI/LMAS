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
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE: {}".format(ASSEMBLY_STATS_GLOBAL_FILE))
    logger.debug("ASSEMBLY_STATS_GLOBAL_FILE_JSON: {}".format(ASSEMBLY_STATS_GLOBAL_FILE_JSON))


def main(assembly_stats_global_file, stats_json):

    # Write CSV file report
    """
    data = {}
    for file in assembly_stats_global_file:
        sample = os.path.basename(file).split('_')[0]

        if sample not in data.keys():
            data[sample] = [file]
        else:
            data[sample].append(file)
    print(data)
    print(data.keys())
    for sample in data.keys():
        print(sample)
        print(data[sample])
        with open(sample + '.csv', "w") as csv_file:
            csv_file.write(','.join(['Assembler', 'Contigs', 'Basepairs', 'Max contig size', 'N50',
                                     'contigs>1000bp (%)', 'Basepairs in contigs>1000bp (%)', 'N50 in contigs>1000bp'])
                           + '\\n')
            for item in data[sample]:
                print(item)
                with open(item, 'r') as stats_file:
                    data = stats_file.read().replace('\\n', '')
                    print(data)
                    csv_file.write(data + '\\n')
    """
    # Write JSON file report
    with open("global_assembly_stats.json", "w") as json_report:
        main_json = {}
        for data_report in stats_json:
            with open(data_report) as f:
                json_data = json.load(f)
                print(json_data)
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
        json_report.write(json.dumps(main_json, separators=(",", ":")))


if __name__ == '__main__':
    #main(ASSEMBLY_STATS_GLOBAL_FILE, ASSEMBLY_STATS_GLOBAL_FILE_JSON)
    main(['mockSample_Pandaseq_global_assembly_stats_global.csv', 'mockSample_MINIA_global_assembly_stats_global.csv', 'mockSample_BCALM2_global_assembly_stats_global.csv', 'ERR2984773_MINIA_global_assembly_stats_global.csv', 'ERR2984773_BCALM2_global_assembly_stats_global.csv', 'mockSample_SKESA_global_assembly_stats_global.csv', 'mockSample_MEGAHIT_global_assembly_stats_global.csv', 'mockSample_IDBA-UD_global_assembly_stats_global.csv', 'ERR2984773_SKESA_global_assembly_stats_global.csv', 'ERR2935805_MINIA_global_assembly_stats_global.csv', 'ERR2984773_Pandaseq_global_assembly_stats_global.csv', 'ERR2984773_MEGAHIT_global_assembly_stats_global.csv', 'ERR2935805_BCALM2_global_assembly_stats_global.csv', 'mockSample_metaSPAdes_global_assembly_stats_global.csv', 'mockSample_VelvetOptimizer_global_assembly_stats_global.csv', 'mockSample_SPAdes_global_assembly_stats_global.csv', 'ERR2935805_SKESA_global_assembly_stats_global.csv', 'ERR2984773_metaSPAdes_global_assembly_stats_global.csv', 'ERR2935805_MEGAHIT_global_assembly_stats_global.csv', 'mockSample_GATBMiniaPipeline_global_assembly_stats_global.csv', 'ERR2984773_GATBMiniaPipeline_global_assembly_stats_global.csv', 'ERR2984773_IDBA-UD_global_assembly_stats_global.csv', 'mockSample_Unicycler_global_assembly_stats_global.csv', 'ERR2935805_IDBA-UD_global_assembly_stats_global.csv', 'ERR2935805_metaSPAdes_global_assembly_stats_global.csv', 'ERR2935805_SPAdes_global_assembly_stats_global.csv', 'ERR2935805_GATBMiniaPipeline_global_assembly_stats_global.csv', 'ERR2984773_SPAdes_global_assembly_stats_global.csv', 'ERR2984773_VelvetOptimizer_global_assembly_stats_global.csv', 'ERR2984773_Unicycler_global_assembly_stats_global.csv', 'ERR2935805_Unicycler_global_assembly_stats_global.csv', 'ERR2935805_VelvetOptimizer_global_assembly_stats_global.csv'], ['mockSample_Pandaseq_report.json', 'mockSample_MINIA_report.json', 'mockSample_BCALM2_report.json', 'ERR2984773_MINIA_report.json', 'ERR2984773_BCALM2_report.json', 'mockSample_SKESA_report.json', 'mockSample_MEGAHIT_report.json', 'mockSample_IDBA-UD_report.json', 'ERR2984773_SKESA_report.json', 'ERR2935805_MINIA_report.json', 'ERR2984773_Pandaseq_report.json', 'ERR2984773_MEGAHIT_report.json', 'ERR2935805_BCALM2_report.json', 'mockSample_metaSPAdes_report.json', 'mockSample_VelvetOptimizer_report.json', 'mockSample_SPAdes_report.json', 'ERR2935805_SKESA_report.json', 'ERR2984773_metaSPAdes_report.json', 'ERR2935805_MEGAHIT_report.json', 'mockSample_GATBMiniaPipeline_report.json', 'ERR2984773_GATBMiniaPipeline_report.json', 'ERR2984773_IDBA-UD_report.json', 'mockSample_Unicycler_report.json', 'ERR2935805_IDBA-UD_report.json', 'ERR2935805_metaSPAdes_report.json', 'ERR2935805_SPAdes_report.json', 'ERR2935805_GATBMiniaPipeline_report.json', 'ERR2984773_SPAdes_report.json', 'ERR2984773_VelvetOptimizer_report.json', 'ERR2984773_Unicycler_report.json', 'ERR2935805_Unicycler_report.json', 'ERR2935805_VelvetOptimizer_report.json'])
