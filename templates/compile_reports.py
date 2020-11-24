#!/usr/bin/env python3

import os
import json
import zipfile
import csv
import re
import fnmatch
from time import gmtime, strftime
try:
    import utils
except ImportError:
    from templates import utils

ASSEMBLY_STATS_REPORT = "$global_assembly_stats"
MAIN_JS = "${js}"
PIPELINE_STATS = "${pipeline_stats}"
CONTIG_SIZE_DISTRIBUTION = "${contig_size_distribution}".split()
MAPPING_STATS_REPORT = "$mapping_assembly_stats"
COMPLETNESS_JSON = "$completness_plots"

ASSEMBLER_PROCESS_LIST = ["BCALM2", "GATBMINIAPIPELINE", "MINIA", "MEGAHIT", "METASPADES", "UNICYCLER", "SPADES",
                          "SKESA", "PANDASEQ", "VELVETOPTIMIZER", "IDBA"]

html_template = """
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>LMAS Report</title>
  </head>
  <body style="background-color: #666666">
    <div id="root"></div>
    <script> const _assemblerPerformanceData = {0} </script>
    <script> const _mainData = {1} </script>
    <script src="./main.js"></script>
  </body>
</html>
"""

#<script> const _fileReportData = {1} </script>

logger = utils.get_logger(__file__)


def _size_coverter(s):
    """Converts size string into megabytes
    Parameters
    ----------
    s : str
        The size string can be '30KB', '20MB' or '1GB'
    Returns
    -------
    float
        With the size in bytes
    """

    if s.upper().endswith("KB"):
        return float(s.rstrip("KB")) / 1024

    elif s.upper().endswith(" B"):
        return float(s.rstrip("B")) / 1024 / 1024

    elif s.upper().endswith("MB"):
        return float(s.rstrip("MB"))

    elif s.upper().endswith("GB"):
        return float(s.rstrip("GB")) * 1024

    elif s.upper().endswith("TB"):
        return float(s.rstrip("TB")) * 1024 * 1024

    else:
        return float(s)


def _hms(s):
    """Converts a hms string into seconds.
    Parameters
    ----------
    s : str
        The hms string can be something like '20s', '1m30s' or '300ms'.
    Returns
    -------
    float
        Time in seconds.
    """

    if s == "-":
        return 0

    if s.endswith("ms"):
        return float(s.rstrip("ms")) / 1000

    fields = list(map(float, re.split("[dhms]", s)[:-1]))
    if len(fields) == 4:
        return fields[0] * 24 * 3600 + fields[1] * 3600 + fields[2] * 60 + \
               fields[3]
    if len(fields) == 3:
        return fields[0] * 3600 + fields[1] * 60 + fields[2]
    elif len(fields) == 2:
        return fields[0] * 60 + fields[1]
    else:
        return fields[0]


def _cpu_load_parser(cpus, cpu_per, t):
    """Parses the cpu load from the number of cpus and its usage
    percentage and returns the cpu/hour measure
    Parameters
    ----------
    cpus : str
        Number of cpus allocated.
    cpu_per : str
        Percentage of cpu load measured (e.g.: 200,5%).
    t : str
        The time string can be something like '20s', '1m30s' or '300ms'.
    """
    try:
        _cpus = float(cpus)
        _cpu_per = float(cpu_per.replace(",", ".").replace("%", ""))
        hours = _hms(t) / 60 / 24

        return ((_cpu_per / (100 * _cpus)) * _cpus) * hours

    except ValueError as e:
        return 0


def _size_compress(s):
    """Shortens a megabytes string.
    """

    if s / 1024 > 1:
        return "{}GB".format(round(s / 1024, 1))
    else:
        return "{}MB".format(s)


def process_performance_data(pipeline_stats):
    """

    :param pipeline_stats:
    :return:
    """
    # Parse performance data
    performance = {}
    with open(pipeline_stats, "r") as pipeline_stats_file:
        csvreader = csv.reader(pipeline_stats_file, delimiter='\t')
        for row in csvreader:
            if row[2] in ASSEMBLER_PROCESS_LIST:
                if row[2] not in performance.keys():
                    performance[row[2]] = {"cpus": [_cpu_load_parser(row[8], row[15], row[13])],
                                           "realtime": [_hms(row[13])],
                                           "rss": [_size_coverter(row[17])],
                                           "rchar": [_size_coverter(row[19])],
                                           "wchar": [_size_coverter(row[20])]}
                else:
                    performance[row[2]]["cpus"].append(_cpu_load_parser(row[8], row[15], row[13]))
                    performance[row[2]]["realtime"].append(_hms(row[13]))
                    performance[row[2]]["rss"].append(_size_coverter(row[17]))
                    performance[row[2]]["rchar"].append(_size_coverter(row[19]))
                    performance[row[2]]["wchar"].append(_size_coverter(row[20]))

    performance_metadata = []

    id_int = 1
    for process_id in performance.keys():
        # time
        time_array = performance[process_id]["realtime"]
        mean_time = round(sum(time_array) / len(time_array), 1)
        mean_time_str = strftime('%H:%M:%S', gmtime(mean_time))

        # cumulative cpu / hours
        cpu_hour = round(sum(performance[process_id]["cpus"]), 2)

        # maximum memory
        max_rss = round(max(performance[process_id]["rss"]))
        rss_str = _size_compress(max_rss)

        # average read size
        avg_rchar = round(sum(performance[process_id]["rchar"]) / len(performance[process_id]["rchar"]))
        rchar_str = _size_compress(avg_rchar)

        # average write size
        avg_wchar = round(sum(performance[process_id]["wchar"]) / len(performance[process_id]["wchar"]))
        wchar_str = _size_compress(avg_wchar)

        performance_metadata.append({"id": id_int, "assembler": process_id, "avgTime": mean_time_str, "cpus": cpu_hour,
                                     "max_rss": rss_str, "avgRead": rchar_str, "avgWrite": wchar_str})
        id_int += 1

    logger.debug("Performance dictionary: {}".format(performance_metadata))

    return performance_metadata


def main(main_js, pipeline_stats, assembly_stats_report, contig_size_plots, mapping_stats_report, completness_plot):

    metadata = {
        "nfMetadata": {
            "scriptId": "${workflow.scriptId}",
            "scriptName": "${workflow.scriptName}",
            "profile": "${workflow.profile}",
            "container": "${workflow.container}",
            "containerEngine": "${workflow.containerEngine}",
            "commandLine": "${workflow.commandLine}",
            "runName": "${workflow.runName}",
            "sessionId": "${workflow.sessionId}",
            "projectDir": "${workflow.projectDir}",
            "launchDir": "${workflow.launchDir}",
            "startTime": "${workflow.start}"
        }
    }

    # Add nextflow metadata
    storage = []
    storage.append(metadata)

    # Assembler performance data:
    performance_metadata = process_performance_data(pipeline_stats)

    # LMAS report

    # main report skeleton
    main_data_js = {}

    # add global stats
    with open(assembly_stats_report) as f:
        assembly_stats_json = json.load(f)
        for sample_id in assembly_stats_json.keys():
            main_data_js[sample_id] = assembly_stats_json[sample_id]

    # add mapping stats
    with open(mapping_stats_report) as f:
        mapping_stats_json = json.load(f)
        for sample_id in mapping_stats_json.keys():
            for reference, mapping_stats_reference in mapping_stats_json[sample_id]["ReferenceTable"].items():
                if reference not in main_data_js[sample_id]["ReferenceTables"].keys():
                    main_data_js[sample_id]["ReferenceTables"][reference] = [mapping_stats_reference]
                else:
                    main_data_js[sample_id]["ReferenceTables"][reference].append(mapping_stats_reference)

    for sample_id in main_data_js.keys():
        # add global plots
        contig_distribution_plot = fnmatch.filter(contig_size_plots, sample_id + '*')[0]
        with open(contig_distribution_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            main_data_js[sample_id]["PlotData"] = {"Global": [plot_json]}

        # add reference plots
        with open(completness_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                reference_plots_json = json.load(reference_plots)
                main_data_js[sample_id]["PlotData"][reference] = reference_plots_json

    logger.debug("Report data dictionary: {}".format(main_data_js))

    with open("pipeline_report.json", "w") as json_fh:
        json_fh.write(json.dumps(main_data_js, separators=(",", ":")))

    with open("pipeline_report.html", "w") as html_fh:
        html_fh.write(html_template.format(performance_metadata, main_data_js))

    with zipfile.ZipFile(main_js) as zf:
        os.mkdir("src")
        zf.extractall("./src")


if __name__ == "__main__":
    main(MAIN_JS, PIPELINE_STATS, ASSEMBLY_STATS_REPORT, CONTIG_SIZE_DISTRIBUTION, MAPPING_STATS_REPORT,
         COMPLETNESS_JSON)
