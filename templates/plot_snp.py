#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to plot the location of snps in the filtered assemblies against the reference sequences.
Produces a line scatter plot for the gap per assembler for each reference.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``snp_coords_dataframes``: list of csv files containing the dataframe with gap location information
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
from collections import Counter
from copy import deepcopy
from itertools import groupby
from plotly.offline import plot
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "17.02.2021"
__template__ = "PLOT_SNP-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    DATAFRAME_LIST = '$snp_coords_dataframes'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DATAFRAME_LIST: {}".format(DATAFRAME_LIST))


def main(dataframes):

    li = []
    for filename in dataframes:
        df = pd.read_csv(filename, index_col=0, header=0)
        li.append(df)

    frame = pd.concat(li, ignore_index=True)

    report_dict = {}
    samples = sorted(frame['Sample'].unique())
    for sample in samples:
        report_dict[sample] = {"PlotData": {}}
        references = sorted(frame['Reference'].unique())
        for reference in references:
            fig = make_subplots(rows=2, cols=1,
                                row_heights=[0.2, 0.8],
                                shared_xaxes=True,
                                vertical_spacing=0.02)

            y = 0
            reference_length = int(
                frame['Reference Length'][frame['Reference'] == reference].unique())

            #indexes = np.arange(reference_length)
            _count = Counter()
            # _count.update({x: 0 for x in indexes}) # initialize all values as 0 counts

            assemblers = sorted(frame['Assembler'].unique(
            ), key=lambda v: v.upper(), reverse=True)
            assemblers_in_plot = []

            for assembler in assemblers:
                coords = frame[(frame['Sample'] == sample) & (frame['Reference'] == reference) &
                               (frame['Assembler'] == assembler)]
                if coords.empty:
                    continue
                else:
                    assemblers_in_plot.append(assembler)
                    coord_list = list(coords['SNP Location'])
                    substitution_list = list(coords['Substitution Type'])

                    # trace with gap location - one per gap
                    fig.add_trace(go.Scattergl(x=coord_list,
                                               y=[y]*len(coord_list),
                                               mode='markers',
                                               marker_symbol='line-ns',
                                               text=substitution_list,
                                               hovertemplate='<b>Substitution:</b>: %{text}' +
                                               '<br><b>Coord</b>: %{x}<br>',
                                               marker=dict(color='#000000', size=12,
                                                           line=dict(width=3,
                                                                     color='#000000')),
                                               name=assembler,
                                               showlegend=False,
                                               ),
                                  row=2, col=1)
                    _count.update(coord_list)
                    y += 1
            # histogram-like plot for snp counts, if snp
            if _count:
                _count_sorted = dict(
                    sorted(_count.items(), key=lambda i: i[0]))
                labels, values = zip(*_count_sorted.items())
                fig.add_trace(go.Bar(x=labels, y=values, marker_color='#000000',
                                     showlegend=False, width=12), row=1, col=1)

            # style plot
            fig.update_xaxes(title_text="{} Bp".format(reference),
                             range=[0, reference_length], row=2, col=1)

            fig.update_yaxes(type='category', tickmode='array',
                             tickvals=list(range(0, y)),
                             ticktext=assemblers_in_plot, row=2, col=1)

            fig.update_layout(plot_bgcolor='rgb(255,255,255)',
                              xaxis=dict(showline=True, zeroline=True,
                                         linewidth=1, linecolor='black',
                                         gridcolor='#DCDCDC'))

            html_filename = '{0}_{1}_snps.html'.format(
                sample, reference.replace(' ', '_'))
            plot(fig, filename=html_filename, auto_open=False)

            plot_json = fig.to_json()

            report_dict[sample]['PlotData'].setdefault(
                reference, []).append(plot_json)

    with open('snps_in_reference.json', 'w') as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(DATAFRAME_LIST)
