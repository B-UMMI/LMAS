#!/usr/bin/env python3
"""
Purpose
-------
Script to obtain information of percentage of mapping contigs/basepairs from the filtered assemblies
and produce a descriptive boxplot for mapped and unmapped contigs per assembler.

Expected input
--------------
The following variables are expected whether using NextFlow or the
:py:func:`main` executor.
- ``coverage_files``: list of files with dataframe as csv to be concatenated and parsed
        e.g.: ``'[SampleA_AssemblerA.csv, SampleB_AssemblerB.csv]'``

Authorship
----------
Inês Mendes, cimendes@medicina.ulisboa.pt
https://github.com/cimendes
"""

import os
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "17.11.2020"
__template__ = "PLOT_CONTIG_DISTRIBUTION-nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    DATAFRAME_FILES = '$dataframes'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("COVERAGE_FILE: {}".format(DATAFRAME_FILES))


def main(dataframe_files):
    """

    :param dataframe_files:
    :return:
    """

    df = pd.concat((pd.read_csv(f) for f in dataframe_files))

    for sample_id in sorted(df['Sample'].unique(), reverse=True):

        fig = go.Figure()

        for assembler in sorted(df['Assembler'].unique(), key=lambda v: v.upper(), reverse=True):

            # mapped contigs as boxplots
            fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] != 'Unmapped') & (df['Assembler'] == assembler) & (df['Sample'] == sample_id)],
                                 name=assembler, boxpoints='outliers',
                                 boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))
            # unmapped contigs as scatter-like plot (boxplot showing only the underlying data)
            fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] == 'Unmapped') & (df['Assembler'] == assembler) & (df['Sample'] == sample_id)],
                                 name=assembler, boxpoints='all', pointpos=0, marker=dict(color='rgba(178,37,34,0.7)'),
                                 line=dict(color='rgba(0,0,0,0)'), fillcolor='rgba(0,0,0,0)'))

        fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Contig size (Log bp)",
                          plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))

        plot(fig, filename='{}_contig_size_distribution.html'.format(sample_id), auto_open=False)
        fig.write_json(file='{}_contig_size_distribution.json'.format(sample_id))


if __name__ == '__main__':
    main(DATAFRAME_FILES)
