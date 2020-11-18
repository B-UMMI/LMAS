#!/usr/bin/python3
import os
import sys
import json
import zipfile
import csv
import re
from datetime import time, timedelta

REPORTS = "${report}".split()
MAIN_JS = "${js}"
PIPELINE_STATS = "${pipeline_stats}"

ASSEMBLER_PROCESS_LIST = ["BCALM2", "GATBMINIAPIPELINE", "MINIA", "MEGAHIT", "METASPADES", "UNICYCLER", "SPADES",
                          "SKESA", "PANDASEQ", "VELVETOPTIMIZER", "IDBA"]

html_template = """
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link href="https://fonts.googleapis.com/icon?family=Material+Icons" rel="stylesheet">
  <title>FlowCraft App</title>
</head>
<body style="background-color: #f2f2f2">
    <div id="app"><!-- React --></div>
</body>
<script> const _assemblerPerformanceData = {0} </script>
<script src="./src/main.js"></script>
</html>
"""

#<script> const _fileReportData = {1} </script>

def detect_time(str):
    """

    :param str:
    :return:
    """
    seconds = re.findall(r"(\d*.\d*)s", str)
    if len(seconds) > 0:
        seconds = int(seconds[0].strip())
    else:
        seconds = 0

    minutes = re.findall(r"(\d*.\d*)m", str)
    if len(minutes) > 0:
        minutes = int(minutes[0].strip())
    else:
        minutes = 0

    hours = re.findall(r"(\d*.\d*)h", str)
    if len(hours) > 0:
        hours = int(hours[0].strip())
    else:
        hours = 0

    return time(hour=hours, minute=minutes, second=seconds)


def conv_MB_to_GB(input_megabyte):
    """

    :param input_megabyte:
    :return:
    """
    gigabyte = float(9.5367431640625E-7)
    convert_gb = gigabyte * input_megabyte
    return convert_gb


def get_max_mem(str):
    """

    :param str:
    :return:
    """
    if 'GB' in str:
        return float(str.replace('GB', '').strip())
    else:
        return conv_MB_to_GB(float(str.replace('MB', '').strip()))


def average(lst):
    return sum(lst) / len(lst)


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


def main(reports, main_js, pipeline_stats):

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

    # Add performance data
    performance = {}
    with open(pipeline_stats, "r") as pipeline_stats_file:
        csvreader = csv.reader(pipeline_stats_file, delimiter='\t')
        for row in csvreader:
            if row[2] in ASSEMBLER_PROCESS_LIST:
                if row[2] not in performance.keys():
                    performance[row[2]] = {"cpus": [_cpu_load_parser(row[8], row[15], row[13])],
                                           "realtime": [detect_time(row[13])],
                                           "rss": [get_max_mem(row[17])],
                                           "rchar": [get_max_mem(row[19])],
                                           "wchar": [get_max_mem(row[20])]}
                else:
                    performance[row[2]]["cpus"].append(_cpu_load_parser(row[8], row[15], row[13]))
                    performance[row[2]]["realtime"].append(detect_time(row[13]))
                    performance[row[2]]["rss"].append(get_max_mem(row[17]))
                    performance[row[2]]["rchar"].append(get_max_mem(row[19]))
                    performance[row[2]]["wchar"].append(get_max_mem(row[20]))

    performance_metadata = []
    id = 1
    for process_id in performance.keys():
        time_list = performance[process_id]["realtime"]
        avg_time = str(timedelta(seconds=sum(map(lambda f: int(f[0])*3600 + int(f[1])*60 + int(f[2]),
                                                 map(lambda f: str(f).split(':'), time_list)))/len(time_list)))
        max_cpus = max(performance[process_id]["cpus"])
        max_rss = max(performance[process_id]["rss"])
        avg_read = average(performance[process_id]["rchar"])
        avg_write = average(performance[process_id]["wchar"])
        performance_metadata.append({"id": id, "assembler": process_id, "avgTime": avg_time, "cpus": max_cpus,
                                     "max_rss": max_rss, "avgRead": avg_read, "avgWrite": avg_write})
        id += 1

    print(performance_metadata)

    for r in reports:
        with open(r) as fh:
            rjson = json.load(fh)
            storage.append(rjson)
            print("{}: {}".format(rjson["processName"],
                                  sys.getsizeof(json.dumps(rjson))))

    with open("pipeline_report.html", "w") as html_fh:
        html_fh.write(html_template.format(performance_metadata))

    """
    with zipfile.ZipFile(main_js) as zf:
        os.mkdir("src")
        zf.extractall("./src")
    """
    """
    with open("pipeline_report.json", "w") as rep_fh:
        rep_fh.write(json.dumps({"data": {"results": storage}},
                                separators=(",", ":")))
    """


if __name__ == "__main__":
    #main(REPORTS, MAIN_JS, PIPELINE_STATS)
    main("", "", "/home/cimendes/Temp/LMAS/pipeline_stats.txt")