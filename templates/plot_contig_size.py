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

    sample_id_list = set()

    for dataframe in dataframe_files:
        sample_name = os.path.basename(dataframe).split('_')[0]
        sample_id_list.add(sample_name)

    df = pd.concat((pd.read_csv(f) for f in dataframe_files))

    print(df)
    sample_id =  list(sample_id_list)[0]

    fig = go.Figure()

    for assembler in sorted(df['Assembler'].unique()):
        contigs = df['Contig Len'][df['Assembler'] == assembler]
        mapped_contigs = df['Contig Len'][(df['Mapped'] != 'Unmapped') & (df['Assembler'] == assembler)]

        print(','.join([assembler, f'{len(mapped_contigs)} ({(len(mapped_contigs) / len(contigs)) * 100:.2f}%)',
                        f'{sum(mapped_contigs)} ({(sum(mapped_contigs) / sum(contigs)) * 100:.2f}%)']))

        # mapped contigs as boxplots
        fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] != 'Unmapped') & (df['Assembler'] == assembler)],
                             name=assembler, boxpoints='outliers',
                             boxmean=False, fillcolor='#D3D3D3', line=dict(color='#000000')))
        # unmapped contigs as scatter-like plot (boxplot showing only the underlying data)
        fig.add_trace(go.Box(x=df['Contig Len'][(df['Mapped'] == 'Unmapped') & (df['Assembler'] == assembler)],
                             name=assembler, boxpoints='all', pointpos=0, marker=dict(color='rgba(178,37,34,0.7)'),
                             line=dict(color='rgba(0,0,0,0)'), fillcolor='rgba(0,0,0,0)'))

    fig.update_layout(showlegend=False, xaxis_type="log", xaxis_title="Contig size (Log bp)",
                      title="Contig size distribution per assembler (contigs over 1000 bp)",
                      plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))

    plot(fig, filename='{}_contig_size_distribution.html'.format(sample_id), auto_open=False)
    fig.write_json(file='{}_contig_size_distribution.json'.format(sample_id))


if __name__ == '__main__':
    main(DATAFRAME_FILES)
