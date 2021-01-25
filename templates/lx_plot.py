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
    LX_FILES = '$lx_files '.split()
    L_TARGET = "$params.l_target"
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("C90_FILES: {}".format(LX_FILES))
    logger.debug("L_TARGET: {}".format(L_TARGET))


def main(c90files, l_target):

    df_Lx = pd.DataFrame(columns=['Reference', 'Assembler', 'Lx', 'nContigs'])

    for file_c90 in c90files:
        sample_name = os.path.basename(file_c90).split('_')[0]
        with open(file_c90) as fh:
            next(fh)  # skip header line
            for line in fh:
                line = line.split(',')
                reference = line[1]
                assembler = line[2]
                Lx = line[3]
                contigs = line[4]
                df_Lx = df_Lx.append({'Sample': sample_name, 'Reference': reference,
                                      'Assembler': assembler, 'Lx': Lx, 'nContigs': contigs}, ignore_index=True)

    # create percentage instead of float
    #df_Lx['Lx'] = df_Lx['Lx'] * 100

    # Create plot - Lx per reference for each sample
    report_dict = {}
    for sample in sorted(df_Lx['Sample'].unique()):
        for reference in sorted(df_Lx['Reference'].unique()):
            fig_Lx = go.Figure()
            i = 0
            for assembler in sorted(df_Lx['Assembler'].unique()):
                fig_Lx.add_trace(go.Scatter(x=df_Lx['Lx'][(df_Lx['Sample'] == sample) &
                                                          (df_Lx['Reference'] == reference) &
                                                          (df_Lx['Assembler'] == assembler)],
                                            y=df_Lx['nContigs'][(df_Lx['Sample'] == sample) &
                                                                (df_Lx['Reference'] == reference) &
                                                                (df_Lx['Assembler'] == assembler)],
                                            name=assembler, line=dict(color=utils.COLOURS[i], width=2)))
                i += 1
                # add target line
            fig_Lx.add_shape(type="line", yref="paper",
                                x0=l_target, y0=0, x1=l_target, y1=1,
                                line=dict(color="#D3D3D3", width=4,dash="dashdot"))

            fig_Lx.update_layout(xaxis_title="L(x) %",
                                 yaxis_title='Contigs',
                                 plot_bgcolor='rgb(255,255,255)',
                                 xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                            gridcolor='#DCDCDC'))

            fig_Lx.update_layout(title={'text': "Lx metric for {}".format(reference),
                                        'y': 1,
                                        'x': 0.5,
                                        'xanchor': 'center',
                                        'yanchor': 'top'})

            fig_Lx.update_layout(updatemenus=list([dict(active=1, buttons=list([
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

            fig_Lx.update_layout(annotations=[dict(text="y axis scale:", x=0, xref="paper", y=1.1, yref="paper",
                                                   align="left", showarrow=False, yanchor="top")])

            plot(fig_Lx, filename='{0}_{1}_lx.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_Lx.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

        print(report_dict[sample]['PlotData'].keys())

    with open("lx.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(LX_FILES, L_TARGET)
