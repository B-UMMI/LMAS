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
    C90_FILES = '$c90_files'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("C90_FILES: {}".format(C90_FILES))


def main(c90files):
    """

    :param c90files:
    :return:
    """

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
            i=0
            for assembler in sorted(df_Lx['Assembler'].unique()):
                fig_Lx.add_trace(go.Scatter(x=df_Lx['Lx'][(df_Lx['Sample'] == sample) &
                                                          (df_Lx['Reference'] == reference) &
                                                          (df_Lx['Assembler'] == assembler)],
                                            y=df_Lx['nContigs'][(df_Lx['Sample'] == sample) &
                                                                (df_Lx['Reference'] == reference) &
                                                                (df_Lx['Assembler'] == assembler)],
                                            name=assembler, line=dict(color=utils.COLOURS[i], width=2)))
                i += 1


            fig_Lx.update_layout(title="Lx metric for {}".format(reference),
                                 xaxis_title="L(x) %",
                                 yaxis_title='Contigs',
                                 plot_bgcolor='rgb(255,255,255)',
                                 xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                            gridcolor='#DCDCDC'))

            plot(fig_Lx, filename='{0}_{1}_c90.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_Lx.to_json()

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
