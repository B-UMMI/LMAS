#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to plot the location of gaps sizes in the filtered assemblies against the reference sequences.
Produces a line scatter plot for the gap per assembler for each reference.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``gap_coords_dataframes``: list of csv files containing the dataframe with gap location information
        e.g.: ``'[SampleA.csv, SampleB.csv]'``

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
import numpy as np
from plotly.offline import plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import Counter
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "15.12.2020"
__template__ = "PLOT_GAP_BOXPLOT-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    DATAFRAME_LIST = '$gap_coords_dataframes'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DATAFRAME_LIST: {}".format(DATAFRAME_LIST))


def main(dataframes):

    li = []

    for filename in dataframes:
        df = pd.read_csv(filename, index_col=0, header=0)
        li.append(df)

    frame = pd.concat(li, ignore_index=True)

    #print(frame)

    report_dict = {}

    for sample in sorted(frame['Sample'].unique()):
        for reference in sorted(frame['Reference'].unique()):
            fig = make_subplots(rows=2, cols=1, row_heights=[0.2, 0.8], shared_xaxes=True, vertical_spacing=0.02)
            y = 0
            indexes = np.arange(frame['Reference Length'][frame['Reference'] == reference].unique())
            width = 1
            _count = Counter()
            _count.update({x: 0 for x in indexes})  # initialize all values as 0 counts
            for assembler in sorted(frame['Assembler'].unique()):
                coords = frame[(frame['Sample'] == sample) & (frame['Reference'] == reference) &
                               (frame['Assembler'] == assembler)]
                for i, row in coords.iterrows():
                    # trace with gap location - one per gap
                    fig.add_trace(go.Scatter(x=[coords.at[i, 'Gap Start'], coords.at[i, 'Gap End']],
                                             y=[y, y], mode='lines', line=dict(color='#000000', width=12),
                                             name=assembler, showlegend=False), row=2, col=1)
                    _count.update(np.arange(coords.at[i, 'Gap Start'], coords.at[i, 'Gap End']+1))
                y += 1

            # histogram-like plot for gap counts
            labels, values = zip(*_count.items())
            fig.add_trace(go.Scatter(x=labels, y=values, mode='lines', line=dict(color='#000000', width=2),
                                     showlegend=False), row=1, col=1)

            # style plot
            fig.update_xaxes(title_text="{} Bp".format(reference), row=2, col=1)
            fig.update_yaxes(type='category', tickmode='array',
                             tickvals=list(range(len(sorted(frame['Assembler'].unique())))),
                             ticktext=sorted(frame['Assembler'].unique()), row=2, col=1)

            fig.update_layout(title="Gaps for {}".format(reference),
                              plot_bgcolor='rgb(255,255,255)',
                              xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                         gridcolor='#DCDCDC'))

            plot(fig, filename='{0}_{1}_gaps.html'.format(sample, reference.replace(' ', '_')), auto_open=False)

            plot_json = fig.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_json]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_json]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_json)

    with open("gaps_in_reference.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(DATAFRAME_LIST)
    #main(["mockSample_IDBA-UD_gaps.csv", "mockSample_MEGAHIT_gaps.csv"])

