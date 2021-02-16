#!/usr/bin/env python3

"""
Purpose
-------
This module is intended parse the results of the misassembly.py for one or more
samples and produce a plot for each.

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
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import pickle
from plotly.subplots import make_subplots
from pandas.core.common import flatten
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "25.01.2021"
__template__ = "PROCESS_MISASSEMBLY-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    MISASSEMBLY_TRACE = '$misassembly_trace'.split()
    MISASSEMBLY_CONTIGS = "$misassembly_contigs".split()
    REPORT_DATA = "$report_data".split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("MISASSEMBLY_TRACE: {}".format(MISASSEMBLY_TRACE))
    logger.debug("MISASSEMBLY_CONTIGS: {}".format(MISASSEMBLY_CONTIGS))
    logger.debug("REPORT_DATA: {}".format(REPORT_DATA))


def main(misassembly_trace, misassembly_contigs, report_data):

    # PLOT
    data_dict = {}
    contig_size = {}

    sorted_misassembly_trace = sorted(misassembly_trace, key=lambda v: v.upper(), reverse=True)

    for file_trace in misassembly_trace:
        sample_name = file_trace.split("_")[0]
        assembler_name = file_trace.split("_")[1]
        with open(file_trace, 'rb') as f:
            trace = pickle.load(f)
        if sample_name not in data_dict.keys():
            data_dict[sample_name] = {assembler_name: trace}
        else:
            if assembler_name not in data_dict[sample_name].keys():
                data_dict[sample_name][assembler_name] = trace
        
    for contig_file in misassembly_contigs:
        sample_name = contig_file.split("_")[0]
        assembler_name = contig_file.split("_")[1]
        with open(contig_file, 'rb') as f:
            x = pickle.load(f)
        if sample_name not in contig_size.keys():
            contig_size[sample_name] = [x]
        else:
            contig_size[sample_name].append(x)


    for sample in data_dict.keys():

        fig = make_subplots(rows=2, cols=1,
                            shared_xaxes=True,
                            row_heights=[0.8, 0.2])
        flatlist = list(flatten(contig_size[sample]))
        fig.add_trace(go.Box(x=flatlist, name="", showlegend=False), row=2, col=1)
        for assembler, trace in sorted(data_dict[sample].items(), key=lambda item: item[0].isupper()):
            logger.debug("Processing {}...".format(assembler))
            fig.add_trace(trace, row=1, col=1)
    
        fig.update_traces(marker=dict(line_width=1, symbol='circle', size=16), col=1)
        fig.update_xaxes(type="log", title="Contig Length")
        fig.update_yaxes(title_text="Number of fragments", row=1, col=1)

        fig.update_layout(title="Misassembled contigs", plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))

        plot(fig, filename='{}_misassembly.html'.format(sample), auto_open=False)
        fig.write_json(file='{}_misassembly.json'.format(sample))

    # TABLE DATA
    master_report_data = {}
    for file_report in report_data:
        with open(file_report) as json_fh:
            data_json = json.load(json_fh)
            print(data_json)
            if data_json["sample"] not in master_report_data.keys():
                master_report_data[data_json["sample"]] = {data_json["assembler"]: data_json["misassembled_contigs"]}
            else:
                master_report_data[data_json["sample"]][data_json["assembler"]] = data_json["misassembled_contigs"]

    with open("misassembly_report.json", "w") as json_report:
        json_report.write(json.dumps(master_report_data, separators=(",", ":")))

if __name__ == '__main__':
    main(MISASSEMBLY_TRACE, MISASSEMBLY_CONTIGS, REPORT_DATA)
