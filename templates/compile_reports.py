#!/usr/bin/env python3

import os
import json
import zipfile
import csv
import re
import fnmatch
from time import gmtime, strftime
from itertools import groupby
from pandas.core.common import flatten
try:
    import utils
except ImportError:
    from templates import utils

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    READS_NUMBER = "$reads_json".split()
    ASSEMBLY_STATS_REPORT = "$global_assembly_stats"
    MAIN_JS = "${js}"
    LMAS_LOGO = "$lmas_png"
    PIPELINE_STATS = "${pipeline_stats}"
    CONTIG_SIZE_DISTRIBUTION = "${contig_size_distribution}".split()
    MAPPING_STATS_REPORT = "$mapping_assembly_stats"
    COMPLETNESS_JSON = "$completness_plots"
    REFERENCE_FILE = "$reference_file"
    LX_JSON = "$lx_plots"
    NAX_JSON = "$nax_plots"
    NGX_JSON = "$ngx_plots"
    SHRIMP_JSON = "$shrimp_plots"
    GAP_REFERENCE_JSON = "$gap_reference_json"
    SNP_REFERENCE_JSON = "$snp_reference_json"
    GAP_HISTOGRAM = "$gap_histogram".split()
    MISASSEMBLY_PLOT = "$plot_misassemblies".split()
    MISASSEMBLY_REPORT = "$misassembly_data"
    MIN_CONTIG_SIZE = "$params.minLength"
    VERSIONS_JSON = "$versions_json"
    MISASSEMBLY_PER_REF = "$misassembly_per_ref"

    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("READS_NUMBER: {}".format(READS_NUMBER))
    logger.debug("ASSEMBLY_STATS_REPORT: {}".format(ASSEMBLY_STATS_REPORT))
    logger.debug("MAIN_JS: {}".format(MAIN_JS))
    logger.debug("MAPPING_STATS_REPORT: {}".format(MAPPING_STATS_REPORT))
    logger.debug("LMAS_LOGO: {}".format(LMAS_LOGO))
    logger.debug("PIPELINE_STATS: {}".format(PIPELINE_STATS))
    logger.debug("CONTIG_SIZE_DISTRIBUTION: {}".format(
        CONTIG_SIZE_DISTRIBUTION))
    logger.debug("COMPLETNESS_JSON: {}".format(COMPLETNESS_JSON))
    logger.debug("REFERENCE_FILE: {}".format(REFERENCE_FILE))
    logger.debug("LX_JSON: {}".format(LX_JSON))
    logger.debug("NAX_JSON: {}".format(NAX_JSON))
    logger.debug("NGX_JSON: {}".format(NGX_JSON))
    logger.debug("SHRIMP_JSON: {}".format(SHRIMP_JSON))
    logger.debug("GAP_REFERENCE_JSON: {}".format(GAP_REFERENCE_JSON))
    logger.debug("SNP_REFERENCE_JSON: {}".format(SNP_REFERENCE_JSON))
    logger.debug("GAP_HISTOGRAM: {}".format(GAP_HISTOGRAM))
    logger.debug("MISASSEMBLY_PLOT: {}".format(MISASSEMBLY_PLOT))
    logger.debug("MISASSEMBLY_REPORT: {}".format(MISASSEMBLY_REPORT))
    logger.debug("MIN_CONTIG_SIZE: {}".format(MIN_CONTIG_SIZE))
    logger.debug("VERSIONS_JSON: {}".format(VERSIONS_JSON))
    logger.debug("MISASSEMBLY_PER_REF: {}".format(MISASSEMBLY_PER_REF))

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
    <script> const _referenceData = {1} </script>
    <script> const _sampleData = {2} </script>
    <script> const _mainDataTables = {3} </script>
    <script> const _mainDataPlots = {4} </script>
    <script> const _sampleList = {5} </script>
    <script> const _minContigSize = {6} </script>
    <script src="./main.js"></script>
  </body>
</html>
"""

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


def process_performance_data(pipeline_stats, versions_json):
    """

    :param pipeline_stats:
    :param versions_json:
    :return:
    """

    # parse assembler versions
    with open(versions_json) as f:
        assembler_versions = json.load(f)

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
                                           "wchar": [_size_coverter(row[20])], }
                else:
                    performance[row[2]]["cpus"].append(
                        _cpu_load_parser(row[8], row[15], row[13]))
                    performance[row[2]]["realtime"].append(_hms(row[13]))
                    performance[row[2]]["rss"].append(_size_coverter(row[17]))
                    performance[row[2]]["rchar"].append(
                        _size_coverter(row[19]))
                    performance[row[2]]["wchar"].append(
                        _size_coverter(row[20]))

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
        avg_rchar = round(
            sum(performance[process_id]["rchar"]) / len(performance[process_id]["rchar"]))
        rchar_str = _size_compress(avg_rchar)

        # average write size
        avg_wchar = round(
            sum(performance[process_id]["wchar"]) / len(performance[process_id]["wchar"]))
        wchar_str = _size_compress(avg_wchar)

        performance_metadata.append({"id": id_int, "assembler": process_id, "version": assembler_versions[process_id],
                                     "avgTime": mean_time_str, "cpus": cpu_hour, "max_rss": rss_str,
                                     "avgRead": rchar_str, "avgWrite": wchar_str})
        id_int += 1

    logger.debug("Performance dictionary: {}".format(performance_metadata))

    return performance_metadata


def process_reference_data(reference_file):
    """
    """

    reference_file_name = os.path.splitext(os.path.basename(reference_file))[0]

    ref_sequences_header = {}

    with open(reference_file) as fh:
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:
            headerStr = utils.REFERENCE_DIC[header.__next__()[1:].strip().split()[
                0]]
            seq = "".join(s.strip() for s in faiter.__next__())
            gc_content = float(
                (seq.count('G') + seq.count('C'))) / len(seq) * 100
            ref_sequences_header[headerStr] = {
                "size": len(seq)/3, "GC": gc_content}

    return_dict = {reference_file_name: ref_sequences_header}

    return return_dict


def process_sample_reads(reads_jsons):
    """

    """

    reads_report = {}

    for sample in reads_jsons:
        with open(sample) as f:
            reads_number_dict = json.load(f)
            for k, v in reads_number_dict.items():
                reads_report[k] = v
    return reads_report


def main(main_js, pipeline_stats, assembly_stats_report, contig_size_plots, mapping_stats_report, completness_plot,
         lmas_logo, reference_file, lx_json, shrimp_json, gap_reference_json, gap_histogram, plot_misassembly, misassembly_report,
         min_contig_size, nax_json, ngx_json, reads_json, snp_reference_json, versions_json, misassembly_per_ref):

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
    logger.debug('Processing {0} data...'.format(pipeline_stats))
    performance_metadata = process_performance_data(
        pipeline_stats, versions_json)

    # Reference info
    reference_info = process_reference_data(reference_file)

    # Sample info
    sample_reads = process_sample_reads(reads_json)

    ################
    # LMAS report

    # main report skeleton
    main_data_tables_js = {}
    main_data_plots_js = {}

    # add global stats
    logger.debug('Processing {0} data...'.format(assembly_stats_report))
    with open(assembly_stats_report) as f:
        assembly_stats_json = json.load(f)
        for sample_id in assembly_stats_json.keys():
            main_data_tables_js[sample_id] = assembly_stats_json[sample_id]

    # add misassembly to global stats
    logger.debug('Processing {0} data...'.format(misassembly_report))
    with open(misassembly_report) as f:
        misassembly_json = json.load(f)
        print(misassembly_json)
        for sample_id in misassembly_json.keys():
            for i in range(0, len(main_data_tables_js[sample_id]["GlobalTable"])):
                assembler = main_data_tables_js[sample_id]["GlobalTable"][i]["assembler"]
                main_data_tables_js[sample_id]["GlobalTable"][i]["filtered"][
                    "misassembled_contigs"] = misassembly_json[sample_id][assembler]

    # add mapping stats
    logger.debug('Processing {0} data...'.format(mapping_stats_report))
    with open(mapping_stats_report) as f:
        mapping_stats_json = json.load(f)
        for sample_id in mapping_stats_json.keys():
            for reference, mapping_stats_reference in mapping_stats_json[sample_id]["ReferenceTable"].items():
                if "ReferenceTables" not in main_data_tables_js[sample_id]:
                    main_data_tables_js[sample_id]["ReferenceTables"] = {}
                if reference not in main_data_tables_js[sample_id]["ReferenceTables"].keys():
                    main_data_tables_js[sample_id]["ReferenceTables"][reference] = [
                        mapping_stats_reference]
                else:
                    main_data_tables_js[sample_id]["ReferenceTables"][reference].append(
                        mapping_stats_reference)
    
    # add misassembly stats
    with open(misassembly_per_ref) as misassembly_fh:
        misassembly_stats = json.load(misassembly_fh)
        for sample_id in main_data_tables_js.keys():
            for reference in  main_data_tables_js[sample_id]["ReferenceTables"].keys():
                for item_row in main_data_tables_js[sample_id]["ReferenceTables"][reference]:
                    print(item_row)
                    assembler = item_row['assembler']
                    if misassembly_stats[sample_id][assembler] == [{}] or reference not in misassembly_stats[sample_id][assembler][0]:
                        item_row['misassembled_contigs'] = 0
                    else:
                         item_row['misassembled_contigs'] = misassembly_stats[sample_id][assembler][0][reference]


    for sample_id in main_data_tables_js.keys():
        main_data_plots_js[sample_id] = {}
        main_data_plots_js[sample_id]["PlotData"] = {}

        # add global plots
        main_data_plots_js[sample_id]["PlotData"]["Global"] = {}

        #   contig size boxplot
        contig_distribution_plot = fnmatch.filter(
            contig_size_plots, sample_id + '*')[0]
        logger.debug('Processing {0} data for {1}...'.format(
            contig_distribution_plot, sample_id))
        with open(contig_distribution_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            main_data_plots_js[sample_id]["PlotData"]["Global"]["contig_size"] = plot_json

        #   gap size boxplot
        gap_distribution_plot = fnmatch.filter(
            gap_histogram, sample_id + '*')[0]
        logger.debug('Processing {0} data for {1}...'.format(
            gap_distribution_plot, sample_id))
        with open(gap_distribution_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            main_data_plots_js[sample_id]["PlotData"]["Global"]["gap_size"] = plot_json

        # misassembly plot
        misassembled_contigs_plot = fnmatch.filter(
            plot_misassembly, sample_id + '*')[0]
        logger.debug('Processing {0} data for {1}...'.format(
            misassembled_contigs_plot, sample_id))
        with open(misassembled_contigs_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            main_data_plots_js[sample_id]["PlotData"]["Global"]["misassembly"] = plot_json

        # add reference plots
        #    completness
        logger.debug('Processing {0} data for {1}...'.format(
            completness_plot, sample_id))
        with open(completness_plot) as plot_fh:
            plot_json = json.load(plot_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    main_data_plots_js[sample_id]["PlotData"][reference] = {
                        "completness": reference_plots_json}

        #    Lx
        logger.debug('Processing {0} data for {1}...'.format(
            lx_json, sample_id))
        with open(lx_json) as lx_fh:
            plot_json = json.load(lx_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "lx": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["lx"] = reference_plots_json
        #   NAx
        logger.debug('Processing {0} data for {1}...'.format(
            nax_json, sample_id))
        with open(nax_json) as nax_fh:
            plot_json = json.load(nax_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "nax": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["nax"] = reference_plots_json

        #   NGx
        logger.debug('Processing {0} data for {1}...'.format(
            ngx_json, sample_id))
        with open(ngx_json) as ngx_fh:
            plot_json = json.load(ngx_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "ngx": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["ngx"] = reference_plots_json

        #    phred-like plot
        logger.debug('Processing {0} data for {1}...'.format(
            shrimp_json, sample_id))
        with open(shrimp_json) as phred_fh:
            plot_json = json.load(phred_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "phred": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["phred"] = reference_plots_json

        #   gap plot
        logger.debug('Processing {0} data for {1}...'.format(
            gap_reference_json, sample_id))
        with open(gap_reference_json) as gap_ref_fh:
            plot_json = json.load(gap_ref_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "gaps": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["gaps"] = reference_plots_json

        # SNP plot
        logger.debug('Processing {0} data for {1}...'.format(
            snp_reference_json, sample_id))
        with open(snp_reference_json) as snp_ref_fh:
            plot_json = json.load(snp_ref_fh)
            for reference, reference_plots in plot_json[sample_id]["PlotData"].items():
                for x in reference_plots:
                    reference_plots_json = json.loads(x)
                    if reference not in main_data_plots_js[sample_id]["PlotData"].keys():
                        main_data_plots_js[sample_id]["PlotData"][reference] = {
                            "snps": reference_plots_json}
                    else:
                        main_data_plots_js[sample_id]["PlotData"][reference]["snps"] = reference_plots_json

    #logger.debug("Report data dictionary: {}".format(main_data_plots_js))

    with open("performance_metadata.json", "w") as json_fh:
        json_fh.write(json.dumps(performance_metadata, separators=(",", ":")))

    with open("reference_metadata.json", "w") as json_fh:
        json_fh.write(json.dumps(reference_info, separators=(",", ":")))

    with open("pipeline_report.json", "w") as json_fh:
        json_fh.write(json.dumps(main_data_tables_js, separators=(",", ":")))
        json_fh.write(json.dumps(main_data_plots_js, separators=(",", ":")))

    with open("index.html", "w") as html_fh:
        html_fh.write(html_template.format(json.dumps(performance_metadata),
                                           json.dumps(reference_info),
                                           json.dumps(sample_reads),
                                           json.dumps(main_data_tables_js),
                                           json.dumps(main_data_plots_js),
                                           list(main_data_tables_js.keys()),
                                           min_contig_size))

    with zipfile.ZipFile(main_js) as zf:
        zf.extractall(".")

    with zipfile.ZipFile(lmas_logo) as zf:
        zf.extractall(".")


if __name__ == "__main__":
    main(MAIN_JS, PIPELINE_STATS, ASSEMBLY_STATS_REPORT, CONTIG_SIZE_DISTRIBUTION, MAPPING_STATS_REPORT,
         COMPLETNESS_JSON, LMAS_LOGO, REFERENCE_FILE, LX_JSON, SHRIMP_JSON, GAP_REFERENCE_JSON, GAP_HISTOGRAM,
         MISASSEMBLY_PLOT, MISASSEMBLY_REPORT, MIN_CONTIG_SIZE, NAX_JSON, NGX_JSON, READS_NUMBER, SNP_REFERENCE_JSON,
         VERSIONS_JSON, MISASSEMBLY_PER_REF)
