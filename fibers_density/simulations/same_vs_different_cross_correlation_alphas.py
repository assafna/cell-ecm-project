import os

import numpy as np
import plotly.graph_objs as go

from fibers_density.simulations import same_vs_different_cross_correlation
from libs.simulations import paths
from plotting import save

# by low connectivity
TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
CELLS_DISTANCE = 7
TIME_LAG = 1
TIME_LAG_INDEX = np.where(np.array(same_vs_different_cross_correlation.TIME_LAGS) == TIME_LAG)[0][0]
ALPHAS = [0, 0.25, 0.5, 0.75, 1]
BETA = 1


def main(_low_connectivity=False):
    _same_alphas_arrays = []
    _different_alphas_arrays = []
    _same_alphas_highest = []
    _different_alphas_highest = []

    for _alpha in ALPHAS:
        print('Alpha:', _alpha)
        _, _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, \
            _different_time_lags_highest = same_vs_different_cross_correlation.compute_fibers_densities(
                _alpha=_alpha, _beta=BETA, _low_connectivity=_low_connectivity)

        _same_alphas_arrays.append(_same_time_lags_arrays[TIME_LAG_INDEX])
        _different_alphas_arrays.append(_different_time_lags_arrays[TIME_LAG_INDEX])
        _same_alphas_highest.append(_same_time_lags_highest)
        _different_alphas_highest.append(_different_time_lags_highest)

    # box plots
    for _name, _arrays in zip(['same', 'different'], [_same_alphas_arrays, _different_alphas_arrays]):
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_alpha,
                    boxpoints=False,
                    line={
                        'width': 1
                    },
                    marker={
                        'size': 10,
                        'color': '#2e82bf'
                    },
                    showlegend=False
                ) for _y, _alpha in zip(_arrays, ALPHAS)
            ],
            layout={
                'xaxis': {
                    'title': 'Alpha',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': ALPHAS
                },
                'yaxis': {
                    'title': _name.capitalize() + ' network correlations',
                    'range': [-1, 1.1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_box_low_con_' + str(_low_connectivity) + '_' + _name
        )

    # bar plot
    for _name, _sums in zip(['same', 'different'], [_same_alphas_highest, _different_alphas_highest]):
        _fig = go.Figure(
            data=go.Bar(
                x=ALPHAS,
                y=[_alpha_sums[TIME_LAG_INDEX] / sum(_alpha_sums) for _alpha_sums in _sums],
                marker={
                    'color': '#2e82bf'
                }
            ),
            layout={
                'xaxis': {
                    'title': 'Alpha',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': ALPHAS
                },
                'yaxis': {
                    'title': 'Highest correlations fraction',
                    'range': [0, 1.1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_bar_low_con_' + str(_low_connectivity) + '_' + _name
        )


if __name__ == '__main__':
    main()
