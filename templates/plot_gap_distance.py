#!/usr/bin/env python3
"""

"""

import os
import json
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
try:
    import utils
except ImportError:
    from templates import utils

__version__ = "0.0.1"
__build__ = "15.12.2020"
__template__ = "PLOT_GAP_HISTOGRAM -nf"

logger = utils.get_logger(__file__)

if __file__.endswith(".command.sh"):
    GAP_JSON = '$gap_distance_json'.split()
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("GAP_JSON: {}".format(GAP_JSON))


# colors for each trace
colours = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c',
           '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']


def main(gap_json):
    """

    :param gap_json:
    :return:
    """

    all_data = {}

    for json_file in gap_json:
        with open(json_file) as jfh:
            data = json.load(jfh)
            for k in data.keys():
                if k not in all_data:
                    all_data[k] = data[k]
                else:
                    for assembler, dist_list in data[k].items():
                        all_data[k][assembler] = dist_list

    print(all_data)

    for sample in all_data.keys():
        fig = go.Figure()
        i = 0
        for k, v in all_data[sample].items():
            fig.add_trace(go.Histogram(x=v, name=k, opacity=0.75, marker_color=colours[i]))
            i += 1

        fig.update_layout(title_text='Gap Distance Distribution', xaxis_type="log",
                          xaxis_title_text='Distance between gaps (log)', yaxis_title_text='Count', barmode='stack',
                          plot_bgcolor='rgb(255,255,255)', xaxis=dict(zeroline=False, gridcolor='#DCDCDC'))

        plot(fig, filename='gap_distance_histogram.html', auto_open=False)
        fig.write_json(file='{}_gap_distance_histogram.json'.format(sample))


if __name__ == '__main__':
    main(GAP_JSON)
