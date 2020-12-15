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
__build__ = "15.12.2020"
__template__ = "PROCESS_SHRIMP-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    PHRED_FILES = '$phred_files'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("PHRED_FILES: {}".format(PHRED_FILES))

# colors for each assembler
colours = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
           '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']

def main(phred_files):
    """

    :param phred_files:
    :return:
    """

    df_phred = pd.DataFrame(columns=['Sample', 'Assembler', 'Reference', 'Contig', 'Contig Length',
                                     'Phred Quality Score'])

    for file_phred in phred_files:
        sample_name = os.path.basename(file_phred).split('_')[0]
        with open(file_phred) as fh:
            next(fh)  # skip header line
            for line in fh:
                line = line.split(',')
                assembler = line[1]
                reference = line[2]
                contig = line[3]
                contig_length = line[4]
                phred_score = line[5]
                df_phred = df_phred.append({'Sample': sample_name, 'Assembler': assembler, 'Reference': reference,
                                            'Contig': contig, 'Contig Length': contig_length,
                                            'Phred Quality Score': phred_score}, ignore_index=True)

    # Create plot - C90 per reference for each sample
    report_dict = {}
    for sample in sorted(df_phred['Sample'].unique()):
        for reference in sorted(df_phred['Reference'].unique()):
            fig_phred = go.Figure()
            i = 0

            for assembler in sorted(df_phred['Assembler'].unique()):
                fig_phred.add_trace(go.Scatter(y=df_phred['Phred Quality Score'][(df_phred['Reference'] == reference) &
                                                                                 (df_phred['Assembler'] == assembler) &
                                                                                 (df_phred['Sample'] == sample)],
                                               x=df_phred['Contig Length'][(df_phred['Reference'] == reference) &
                                                                           (df_phred['Assembler'] == assembler) &
                                                                           (df_phred['Sample'] == sample)],
                                    name=assembler,
                                    opacity=0.7,
                                    mode='markers',
                                    marker=dict(color=colours[i], size=12, line=dict(width=1, color='black'))))
                i += 1

            fig_phred.update_layout(title="Phred-like score metric for {}".format(reference),
                                    xaxis_title="Contig size",
                                    yaxis_title="Score",
                                    plot_bgcolor='rgb(255,255,255)',
                                    xaxis=dict(showline=True, zeroline=False, linewidth=1, linecolor='black',
                                               gridcolor='#DCDCDC'))

            plot(fig_phred, filename='{0}_{1}_phred.html'.format(sample, reference.replace(' ', '_')), auto_open=False)

            plot_species = fig_phred.to_json()

            if sample not in report_dict.keys():
                report_dict[sample] = {"PlotData": {reference: [plot_species]}}
            else:
                if reference not in report_dict[sample]["PlotData"].keys():
                    report_dict[sample]["PlotData"][reference] = [plot_species]
                else:
                    report_dict[sample]["PlotData"][reference].append(plot_species)

        print(report_dict[sample]['PlotData'].keys())

    with open("phred.json", "w") as json_report:
        json_report.write(json.dumps(report_dict, separators=(",", ":")))


if __name__ == '__main__':
    main(PHRED_FILES)
