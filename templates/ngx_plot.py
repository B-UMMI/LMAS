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
    NGX_FILES = '$ngx_files '.split()
    N_TARGET = float("$params.n_target")
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("NAX_FILES: {}".format(NGX_FILES))
    logger.debug("N_TARGET: {}".format(N_TARGET))


def main(ngx_files, n_target):

    df_ngx = pd.DataFrame(columns=['Sample', 'Reference', 'Assembler', 'NGx', 'Basepairs'])


    for file_ngx in ngx_files:
        sample_name = os.path.basename(file_ngx).split('_')[0]
        data = pd.read_csv(file_ngx)
        data['Sample'] = sample_name
        df_ngx = pd.concat([df_ngx, data], ignore_index=True)
    
    # Create plot - Lx per reference for each sample
    report_dict = {}
    for sample in sorted(df_ngx['Sample'].unique()):
        for reference in sorted(df_ngx['Reference'].unique()):
            fig_ngx = go.Figure()
            i = 0
            for assembler in sorted(df_ngx['Assembler'].unique(), key=lambda v: v.upper()):
                if df_ngx['Basepairs'][(df_ngx['Sample'] == sample) & (df_ngx['Reference'] == reference) & (df_ngx['Assembler'] == assembler)].nunique() > 1:
                    fig_ngx.add_trace(go.Scatter(x=df_ngx['NGx'][(df_ngx['Sample'] == sample) &
                                                                    (df_ngx['Reference'] == reference) &
                                                                    (df_ngx['Assembler'] == assembler)],
                                                    y=df_ngx['Basepairs'][(df_ngx['Sample'] == sample) &
                                                                        (df_ngx['Reference'] == reference) &
                                                                        (df_ngx['Assembler'] == assembler)],
                                                    name=assembler, line=dict(color=utils.COLOURS[i], width=2)))
                    i += 1
            
            fig_ngx.add_shape(type="line", yref="paper",
                              x0=n_target, y0=0, x1=n_target, y1=1,
                              line=dict(color="#D3D3D3", width=4,dash="dashdot"))

            fig_ngx.update_layout(xaxis_title="NG(x) %",
                                  yaxis_title='Basepairs',
                                  plot_bgcolor='rgb(255,255,255)',
                                  xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                             gridcolor='#DCDCDC'))


            fig_ngx.update_layout(title={
                'text': "NGx metric for {}".format(reference),
                'y': 1,
                'x': 0.5,
                'xanchor': 'center',
                'yanchor': 'top'})

            plot(fig_ngx, filename='{0}_{1}_ngx.html'.format(sample, reference.replace(' ', '_')), auto_open=False)
            plot_species = fig_ngx.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

    with open("ngx.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(NGX_FILES, N_TARGET)

