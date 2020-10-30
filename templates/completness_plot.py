#!/usr/bin/env python3
"""
Purpose
-------
This script takes csv tables by species and plots them as series of scatter plots

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``coverage_files``: list of files with completeness information to be plotted
        e.g.: ``'[SampleA_AssemblerA.csv, SampleB_AssemblerB.csv]'``

Authorship
----------
Rafael Mamede, rmamede@medicina.ulisboa.ptt
https://github.com/rfm-targa
"""

import os
import numpy as np
import pandas as pd
from plotly import subplots
from scipy import interpolate
import plotly.graph_objs as go
from plotly.offline import plot
import utils


__version__ = "0.0.1"
__build__ = "30.10.2020"
__template__ = "PROCESS_COMPLETNESS-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    COVERAGE_FILES = '$coverage_files'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("COVERAGE_FILES: {}".format(COVERAGE_FILES))


def plot_data(sample_id, species_data):
    """

    :param sample_id:
    :param species_data: dictionary
    :return:
    """
    interpolation_xvalues = [0, 40, 80, 160, 320, 640, 1280, 2560]
    interpolation_function = interpolate.interp1d(interpolation_xvalues, np.arange(len(interpolation_xvalues)))

    # colors for each assembler
    colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
              '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
    # create a tracer for each assembler data point
    # simpler to manage colors and legends
    tracers = {}
    i = 0
    j = 0
    group = 1
    for k, v in species_data.items():
        for a, d in v.items():
            legend = True if j == 0 else False
            text = str(d[1]) + '<br>' + a
            tracer = go.Scatter(x=list(interpolation_function([d[1]])),
                                y=[d[0]],
                                name=a,
                                showlegend=legend,
                                legendgroup='group{0}'.format(group),  # group legends
                                opacity=1,
                                mode='markers',
                                marker=dict(color=colors[i],
                                            size=7,
                                            line=dict(width=1, color='black')),
                                text=text,
                                hoverinfo='y+text'
                                )
            i += 1
            group += 1
            tracers.setdefault(k, []).append(tracer)

        i = 0
        j += 1
        group = 1

    num_cols = 2
    num_rows = len(tracers) / num_cols
    num_rows = int(num_rows) if len(tracers) % num_cols == 0 else int(num_rows) + 1
    # define number and organization of subplots
    meta_subplots = subplots.make_subplots(rows=num_rows,
                                           cols=num_cols,
                                           # shared_xaxes=True,
                                           shared_yaxes=True,
                                           subplot_titles=list(species_data.keys()),
                                           horizontal_spacing=0.04,
                                           vertical_spacing=0.1)

    # add assemblers tracers to subplot
    r = 1
    c = 1
    for k, v in tracers.items():
        for t in v:
            meta_subplots.add_trace(t, r, c)
        c += 1
        if c > num_cols:
            r += 1
            c = 1

    # define xaxes attributes
    xaxis_range = list(interpolation_function(interpolation_xvalues))
    meta_subplots.update_xaxes(showgrid=False,
                               showline=True,
                               linecolor='black',
                               linewidth=2,
                               ticks='outside',
                               tickcolor='black',
                               tickwidth=2,
                               tickfont=dict(color='black',
                                             size=12),
                               range=[xaxis_range[0], xaxis_range[-1]],
                               tickvals=xaxis_range,
                               ticktext=interpolation_xvalues)

    # define yaxes attributes
    meta_subplots.update_yaxes(showgrid=False,
                               showline=True,
                               linecolor='black',
                               linewidth=2,
                               tickfont=dict(color='black',
                                             size=12),
                               range=[0, 1.05],
                               tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1])

    # add xaxis title
    xaxis_title = dict(x=0.5, y=-0.1, xref='paper', yref='paper', text='Number of Contigs', showarrow=False)
    yaxis_title = dict(x=-0.07, y=0.5, xref='paper', yref='paper', text='Breadth of Coverage (%)', textangle=-90,
                       showarrow=False)
    meta_subplots.layout.annotations = meta_subplots.layout.annotations + (xaxis_title, yaxis_title)

    # change background color for all subplots
    meta_subplots['layout']['plot_bgcolor'] = 'rgb(255,255,255)'

    # change legend attributes
    meta_subplots['layout']['legend']['font'] = dict(color='black', size=12)
    meta_subplots['layout']['legend']['orientation'] = 'h'
    meta_subplots['layout']['legend']['x'] = 0
    meta_subplots['layout']['legend']['y'] = 1.15
    meta_subplots['layout']['legend']['borderwidth'] = 2

    # change annotations attributes
    for i in meta_subplots['layout']['annotations']:
        i['font']['color'] = 'black'
        i['font']['size'] = 16

    plot(meta_subplots, filename='{}_breadth_of_coverage_plot.html'.format(sample_id), auto_open=False)


def main(coverage_files):

    all_data={}
    for file in coverage_files:
        species_data = {}
        sample_name = os.path.basename(file).split('_')[0]
        assembler_name = os.path.basename(file).split('_')[1]
        logger.debug('Processing {0} {1} data...'.format(sample_name, assembler_name))

        # import table with data
        data = pd.read_csv(file)
        species = list(data['Reference'])
        coverage = list(data['Breadth of Coverage'])
        contigs = list(data['Contigs'])

        for i, s in enumerate(species):
            species_data.setdefault(s, {})[assembler_name] = (coverage[i], contigs[i])

        all_data[sample_name] = species_data

    for sample_id in all_data.keys():
        plot_data(sample_id, all_data[sample_id])


if __name__ == '__main__':
    main(COVERAGE_FILES)