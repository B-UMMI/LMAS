#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to plot the distribution of gap sizes of the filtered assemblies respective to the reference sequences.
Produces a boxplot for the gap size per assembler.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``gap_distance_json``: list of JSON files containing the sample ID and the list of gap sizes
        e.g.: ``'[SampleA.json, SampleB.json]'``

Expected input
--------------
This script takes the following arguments (in this order):
  * List of Paths with the JSON containing the sample ID and the list of gap sizes

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import json
import pandas as pd
from pandas.core.common import flatten
from plotly.offline import plot
import plotly.graph_objects as go
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "15.12.2020"
__template__ = "PLOT_GAP_BOXPLOT-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    GAP_JSON = '$gap_distance_json'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("GAP_JSON: {}".format(GAP_JSON))

COLUMNS = ['Assembler', 'Gap size']  # columns for dataframe


def main(gap_json):

    all_data = {}

    for json_file in gap_json:
        with open(json_file) as jfh:
            data = json.load(jfh)
            for sample in data.keys():
                for assembler in data[sample].keys():
                    if sample not in all_data.keys():
                        all_data[sample] = {assembler: [gap for gap in data[sample][assembler]]}
                        print(all_data[sample][assembler])
                    else:
                        if assembler not in all_data[sample].keys():
                            all_data[sample][assembler] = [gap for gap in data[sample][assembler]]
                        else:
                            all_data[sample][assembler].append(gap for gap in data[sample][assembler])

    for sample in all_data.keys():

        df = pd.DataFrame(columns=COLUMNS)

        fig = go.Figure()
        for k, v in all_data[sample].items():
            flatlist = list(flatten(v))
            for gap in flatlist:
                df = df.append({'Assembler': k, 'Gap size': gap}, ignore_index=True)
        
        for assembler in sorted(df['Assembler'].unique(), key=lambda v: v.upper(), reverse=True):
            fig.add_trace(go.Box(x=df['Gap size'][df['Assembler'] == assembler],
                                 name=assembler, boxpoints='outliers',
                                 boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))

        fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Gap size (Log bp)",
                          title="Gap size distribution per assembler (contigs over 1000 bp)",
                          plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))

        plot(fig, filename='{}_gap_size_boxplot.html'.format(sample), auto_open=False)
        fig.write_json(file='{}_gap_distance_histogram.json'.format(sample))


if __name__ == '__main__':
    main(GAP_JSON)
