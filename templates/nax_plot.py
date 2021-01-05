#!/usr/bin/env python3

"""
"""

import os
import pandas as pd
import json
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
    NAX_FILES = '$nax_files '.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("NAX_FILES: {}".format(NAX_FILES))


def main(nax_files):
    """

    :param nax_files:
    :return:
    """

    df_nax = pd.DataFrame(columns=['Reference', 'Assembler', 'NAx', 'Basepairs'])

    for file_nax in nax_files:
        sample_name = os.path.basename(file_nax).split('_')[0]
        with open(file_nax) as fh:
            next(fh)  # skip header line
            for line in fh:
                line = line.split(',')
                reference = line[1]
                assembler = line[2]
                nax = line[3]
                basepairs = line[4]
                df_nax = df_nax.append({'Sample': sample_name, 'Reference': reference,
                                        'Assembler': assembler, 'NAx': nax, 'Basepairs': basepairs}, ignore_index=True)

    # create percentage instead of float
    #df_Lx['Lx'] = df_Lx['Lx'] * 100

    # Create plot - Lx per reference for each sample
    report_dict = {}
    for sample in sorted(df_nax['Sample'].unique()):
        for reference in sorted(df_nax['Reference'].unique()):
            fig_nax = go.Figure()
            i = 0
            for assembler in sorted(df_nax['Assembler'].unique()):
                fig_nax.add_trace(go.Scatter(x=df_nax['NAx'][(df_nax['Sample'] == sample) &
                                                          (df_nax['Reference'] == reference) &
                                                          (df_nax['Assembler'] == assembler)],
                                             y=df_nax['Basepairs'][(df_nax['Sample'] == sample) &
                                                                (df_nax['Reference'] == reference) &
                                                                (df_nax['Assembler'] == assembler)],
                                             name=assembler, line=dict(color=utils.COLOURS[i], width=2)))
                i += 1

            fig_nax.update_layout(xaxis_title="NA(x) %",
                                  yaxis_title='Basepairs',
                                  plot_bgcolor='rgb(255,255,255)',
                                  xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                             gridcolor='#DCDCDC'))

            fig_nax.update_layout(title={
                'text': "NAx metric for {}".format(reference),
                'y': 1.2,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'})

            fig_nax.update_layout(updatemenus=list([dict(active=1, buttons=list([
                dict(label='Log Scale', method='update',
                     args=[{'visible': [True, True]}, {'yaxis': {'type': 'log'}}]),
                dict(label='Linear Scale', method='update',
                     args=[{'visible': [True, False]}, {'yaxis': {'type': 'linear'}}])
            ]),
                                                         type="buttons",
                                                         direction="right",
                                                         pad={"r": 10, "t": 10},
                                                         showactive=True,
                                                         xanchor="left",
                                                         x=0.05,
                                                         y=1.12,
                                                         yanchor="top")]))

            fig_nax.update_layout(annotations=[dict(text="y axis scale:", x=0, xref="paper", y=1.12, yref="paper",
                                                    align="left", showarrow=False)])

            """
            fig_nax.update_layout(updatemenus=list([dict(active=1, buttons=list([
                dict(label='Log Scale', method='update',
                     args=[{'visible': [True, True]},
                           {'yaxis': {'type': 'log'}}]),
                dict(label='Linear Scale',
                     method='update',
                     args=[{'visible': [True, False]},
                           {'yaxis': {'type': 'linear'}}])
            ]),
                                                         )
                                                    ]))
            """

            plot(fig_nax, filename='{0}_{1}_nax.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_nax.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

        print(report_dict[sample]['PlotData'].keys())

    with open("nax.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(NAX_FILES)
