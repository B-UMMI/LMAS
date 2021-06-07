#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Purpose
-------
Script to plot the location of misassembled contigs in the filtered assemblies against the reference sequences.
Produces a line scatter plot for the misassembled blocks per assembler for each reference.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``misassembly_dataframes``: list of csv files containing the dataframe with gap location information
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
__build__ = "07.06.2021"
__template__ = "PLOT_MISASSEMBLY-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    DATAFRAME_LIST = '$misassembly_dataframes'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DATAFRAME_LIST: {}".format(DATAFRAME_LIST))


def determine_missing_intervals(intervals, total_len):
    """
    """

    start = 0
    missing_regions = []
    for i in intervals:
        diff = i[0] - start
        if diff > 0:
            missing_regions.append([start, start+diff-1])
            start = i[1]

    # add terminal region
    if start != total_len:
        missing_regions.append([start, total_len])

    return missing_regions


def merge_intervals(intervals):
    """ Merges intersecting intervals.
    """
    merged = [deepcopy(intervals[0])]
    for current in intervals[1:]:
        previous = merged[-1]
        # current and previous intervals intersect
        if current[0] <= previous[1]:
            # determine top position
            previous[1] = max(previous[1], current[1])
            # merge coverage dictionaries
            previous_cov = previous[2]
            current_cov = current[2]
            for k, v in current_cov.items():
                if k not in previous_cov:
                    previous_cov[k] = v
                else:
                    previous_cov[k] += v
            previous[2] = previous_cov
        # current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged


def intervals_subgroups(intervals):
    """
    """

    subgroups = {}
    for i in intervals:
        start = i[0]
        current_interval = i[2]
        # identify groups of subsequent equal values
        values_groups = [list(v)
                         for k, v in groupby(current_interval.values())]
        for g in values_groups:
            # keep only the start and end points for each group
            subgroups[start] = g[0]
            subgroups[start+len(g)-1] = g[0]
            start += len(g)

    return subgroups


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
            gaps_intervals = []
            assemblers = sorted(frame['Assembler'].unique(
            ), key=lambda v: v.upper(), reverse=True)
            assemblers_in_plot = []
            for assembler in assemblers:
                print(assembler)
                coords = frame[(frame['Sample'] == sample) & (frame['Reference'] == reference) &
                               (frame['Assembler'] == assembler)]
                if coords.empty:
                    continue
                else:
                    assemblers_in_plot.append(assembler)
                    starts = list(coords['Ref Start'])
                    stops = list(coords['Ref End'])
                    for i in range(len(starts)):
                        gap_size = stops[i]-starts[i]
                        text='<b>Block size</b>: ' + gap_size + '<br> <b>Contig:</b> ' + coords['Contig']
                        # trace with misassembly location - one per gap
                        fig.add_trace(go.Scatter(x=[starts[i], stops[i]],
                                                 y=[y, y],
                                                 mode='lines',
                                                 line=dict(
                                                     color='#000000', width=12),
                                                 name=assembler,
                                                 showlegend=False,
                                                 text=[text, text],
                                                 hovertemplate='%{text}' +
                                                 '<br><b>Coord</b>: %{x}<br>'),
                                      row=2, col=1)
                        gaps_dict = {i: 1 for i in range(
                            starts[i], stops[i]+1)}
                        gaps_intervals.append(
                            [starts[i], stops[i]+1, gaps_dict])
                    y += 1

            reference_length = int(
                frame['Reference Length'][frame['Reference'] == reference].unique())
            if len(gaps_intervals) == 0:
                data_points = gaps_intervals
            else:
                # sort intervals before merging
                gaps_intervals = sorted(gaps_intervals, key=lambda x: x[0])
                merged_intervals = merge_intervals(gaps_intervals)

                # determine missing intervals
                missing_intervals = determine_missing_intervals(
                    merged_intervals, reference_length)

                # identify start and end points for gaps subgroups
                gaps_points = intervals_subgroups(merged_intervals)

                # add points for intervals without gaps
                for i in missing_intervals:
                    gaps_points[i[0]] = 0
                    gaps_points[i[1]] = 0

                data_points = sorted(gaps_points.items(), key=lambda x: x[0])

            labels = [c[0] for c in data_points]
            values = [c[1] for c in data_points]
            # histogram-like plot for gap counts
            fig.add_trace(go.Scatter(x=labels, y=values,
                                     mode='lines',
                                     line=dict(color='#000000', width=0.5),
                                     fill='tozeroy',
                                     showlegend=False),
                          row=1, col=1)

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

            html_filename = '{0}_{1}_misassembly.html'.format(
                sample, reference.replace(' ', '_'))
            plot(fig, filename=html_filename, auto_open=False)

            plot_json = fig.to_json()

            report_dict[sample]['PlotData'].setdefault(
                reference, []).append(plot_json)

    with open('misassembly_in_reference.json', 'w') as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(DATAFRAME_LIST)
