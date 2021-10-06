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
    N_TARGET = float("$params.n_target")
    SCALE = '$scale'
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("NAX_FILES: {}".format(NAX_FILES))
    logger.debug("N_TARGET: {}".format(N_TARGET))
    logger.debug("SCALE: {}".format(SCALE))


def main(nax_files, n_target, scale):

    df_nax = pd.DataFrame(columns=['Sample', 'Reference', 'Assembler', 'NAx', 'Basepairs'])

    for file_nax in nax_files:
        sample_name = os.path.basename(file_nax).split('_')[0]
        data = pd.read_csv(file_nax)
        data['Sample'] = sample_name
        df_nax = pd.concat([df_nax, data], ignore_index=True)

    # Create plot - Lx per reference for each sample
    report_dict = {}
    for sample in sorted(df_nax['Sample'].unique()):
        for reference in sorted(df_nax['Reference'].unique()):
            fig_nax = go.Figure()
            i = 0
            for assembler in sorted(df_nax['Assembler'].unique(), key=lambda v: v.upper()):
                if df_nax['Basepairs'][(df_nax['Sample'] == sample) &(df_nax['Reference'] == reference) & (df_nax['Assembler'] == assembler)].nunique() > 1:
                    fig_nax.add_trace(go.Scatter(x=df_nax['NAx'][(df_nax['Sample'] == sample) &
                                                                (df_nax['Reference'] == reference) &
                                                                (df_nax['Assembler'] == assembler)],
                                                    y=df_nax['Basepairs'][(df_nax['Sample'] == sample) &
                                                                    (df_nax['Reference'] == reference) &
                                                                    (df_nax['Assembler'] == assembler)],
                                                    name=assembler, line=dict(color=utils.COLOURS[i], width=2)))
                    i += 1
            
            fig_nax.add_shape(type="line", yref="paper",
                              x0=n_target*100, y0=0, x1=n_target*100, y1=1,
                              line=dict(color="#D3D3D3", width=4,dash="dashdot"))

            fig_nax.update_layout(xaxis_title="NA(x) %",
                                  yaxis_title='Basepairs',
                                  plot_bgcolor='rgb(255,255,255)',
                                  xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                             gridcolor='#DCDCDC'))
            
            fig_nax.update_xaxes(type=scale)

            plot(fig_nax, filename='{0}_{1}_nax.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_nax.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

    with open("nax.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(NAX_FILES, N_TARGET, SCALE)
