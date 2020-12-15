#!/usr/bin/env python3

"""
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
__build__ = "14.12.2020"
__template__ = "PROCESS_C90-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    C90_FILES = '$c90_files'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("C90_FILES: {}".format(C90_FILES))

# colors for each trace
colours = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
           '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']


def main(c90files):
    """

    :param c90files:
    :return:
    """

    df_c90 = pd.DataFrame(columns=['Sample', 'Reference', 'Assembler', 'C90'])

    for file_c90 in c90files:
        sample_name = os.path.basename(file_c90).split('_')[0]
        with open(file_c90) as fh:
            next(fh)  # skip header line
            for line in fh:
                line = line.split(',')
                reference = line[1]
                assembler = line[2]
                c90 = line[3]
                df_c90 = df_c90.append({'Sample': sample_name, 'Reference': reference,
                                        'Assembler': assembler, 'C90': c90}, ignore_index=True)

    # Create plot - C90 per reference for each sample
    report_dict = {}
    for sample in sorted(df_c90['Sample'].unique()):
        for reference in sorted(df_c90['Reference'].unique()):
            data = [go.Bar(y=df_c90['Assembler'][(df_c90['Sample'] == sample) &
                                                           (df_c90['Reference'] == reference)],
                           x=df_c90['C90'][(df_c90['Sample'] == sample) & (df_c90['Reference'] == reference)],
                           marker_color='lightslategray', orientation='h')]

            updatemenus = list([
                dict(active=1,
                     buttons=list([
                         dict(label='Log Scale',
                              method='update',
                              args=[{'visible': [True, True]},
                                    {'title': 'Log scale',
                                     'xaxis': {'type': 'log'}}]),
                         dict(label='Linear Scale',
                              method='update',
                              args=[{'visible': [True, False]},
                                    {'title': 'Linear scale',
                                     'xaxis': {'type': 'linear'}}])
                     ]),)])

            layout = dict(updatemenus=updatemenus, title='Linear scale')
            fig_c90 = go.Figure(data=data, layout=layout)

            fig_c90.update_layout(title="C90 metric for {}".format(reference),
                                  xaxis_title="Contigs",
                                  plot_bgcolor='rgb(255,255,255)',
                                  xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                             gridcolor='#DCDCDC'))

            plot(fig_c90, filename='{0}_{1}_c90.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_c90.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

        print(report_dict[sample]['PlotData'].keys())

    with open("c90.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(C90_FILES)
