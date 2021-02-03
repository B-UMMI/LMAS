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
Rafael Mamede, rmamede@medicina.ulisboa.pt
https://github.com/rfm-targa
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import numpy as np
import pandas as pd
import json
from scipy import interpolate
import plotly.graph_objs as go
from plotly.offline import plot
try:
    import utils
except ImportError:
    from templates import utils


__version__ = "0.0.1"
__build__ = "30.10.2020"
__template__ = "PROCESS_COMPLETNESS-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    COVERAGE_FILES = '$coverage_files'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("COVERAGE_FILE: {}".format(COVERAGE_FILES))


def plot_data(species_data, sample_id):
    """

    :param sample_id: string
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

    report_json = {}

    for species in species_data.keys():
        i = 0
        to_plot = go.Figure()
        for a, d in species_data[species].items():
            text = str(d[1]) + '<br>' + a
            to_plot.add_trace(go.Scatter(x=list(interpolation_function([d[1]])),
                                         y=[d[0]],
                                         name=a,
                                         showlegend=True,
                                         opacity=1,
                                         mode='markers',
                                         marker=dict(color=colors[i],
                                                     size=12,
                                                     line=dict(width=1, color='black')),
                                         text=text,
                                         hoverinfo='y+text'
                                         ))
            i += 1

        #add title
        to_plot.update_layout(title_text="<br>Needs beteer title here for.. {} ".format(species), title_x=0.5,
                              xaxis_title="Number of contigs", yaxis_title="Breadth of coverage")

        # define xaxes attributes
        xaxis_range = list(interpolation_function(interpolation_xvalues))
        to_plot.update_xaxes(showgrid=False,
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
        to_plot.update_yaxes(showgrid=False,
                             showline=True,
                             linecolor='black',
                             linewidth=2,
                             tickfont=dict(color='black',
                                           size=12),
                             range=[0, 1.05],
                             tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1])

        # change background color for all subplots
        to_plot['layout']['plot_bgcolor'] = 'rgb(255,255,255)'

        # change legend attributes
        to_plot['layout']['legend']['font'] = dict(color='black', size=12)

        # change annotations attributes
        for i in to_plot['layout']['annotations']:
            i['font']['color'] = 'black'
            i['font']['size'] = 16

        plot(to_plot, filename='{0}_{1}_breadth_of_coverage_plot.html'.format(sample_id, species.replace(' ', '_')),
             auto_open=False)

        if sample_id not in report_json.keys():
            report_json[sample_id] = {species: to_plot.to_json()}
        else:
            report_json[sample_id][species] = to_plot.to_json()

    return report_json


def main(coverage_files):

    all_data = {}
    species_data = {}
    for coverage_file in coverage_files:
        sample_name = os.path.basename(coverage_file).split('_')[0]
        assembler_name = os.path.basename(coverage_file).split('_')[1]
        logger.debug('Processing {0} {1} data...'.format(sample_name, assembler_name))

        # import table with data
        data = pd.read_csv(coverage_file)
        species = list(data['Reference'])
        coverage = list(data['Breadth of Coverage'])
        contigs = list(data['Contigs'])

        for i, s in enumerate(species):
            if i not in species_data.keys():
                species_data.setdefault(s, {})[assembler_name] = (coverage[i], contigs[i])

    all_data[sample_name] = species_data
    print(all_data)

    with open("{}_completness_plots.json".format(sample_id), "w") as json_report:
        report_dict = {}
        for sample_id in all_data.keys():
            report_json = plot_data(all_data[sample_id], sample_id)

            for species, plot_species in report_json[sample_id].items():
                if sample_id not in report_dict.keys():
                    report_dict[sample_id] = {"PlotData": {species: [plot_species]}}
                else:
                    if species not in report_dict[sample_id]["PlotData"].keys():
                        report_dict[sample_id]["PlotData"][species] = [plot_species]
                    else:
                        report_dict[sample_id]["PlotData"][species].append(plot_species)

        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(COVERAGE_FILES)

