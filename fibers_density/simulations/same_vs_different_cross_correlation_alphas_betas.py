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
ALPHA = 1
BETAS = [1, 1.05, 1.1, 1.2]
BETA = 1


def main(_type='alpha', _low_connectivity=False, _plots=None):
    _same_arrays = []
    _different_arrays = []
    _same_highest = []
    _different_highest = []

    if _type == 'alpha':
        _alphas = ALPHAS
        _betas = [BETA] * len(ALPHAS)
        _names = _alphas
    elif _type == 'beta':
        _alphas = [ALPHA] * len(BETAS)
        _betas = BETAS
        _names = _betas
    else:
        raise Exception('No such type. Only \'alpha\' or \'beta\' are acceptable.')

    for _alpha, _beta in zip(_alphas, _betas):
        print('Alpha:', _alpha, 'beta:', _beta)
        _, _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, \
            _different_time_lags_highest = same_vs_different_cross_correlation.compute_fibers_densities(
                _alpha=_alpha, _beta=_beta, _low_connectivity=_low_connectivity)

        _same_arrays.append(_same_time_lags_arrays[TIME_LAG_INDEX])
        _different_arrays.append(_different_time_lags_arrays[TIME_LAG_INDEX])
        _same_highest.append(_same_time_lags_highest)
        _different_highest.append(_different_time_lags_highest)

    if _plots is not None:

        # box plot
        for _name, _arrays in zip(['same', 'different'], [_same_arrays, _different_arrays]):
            if _name in _plots:
                _fig = go.Figure(
                    data=[
                        go.Box(
                            y=_y,
                            name=_name,
                            boxpoints=False,
                            line={
                                'width': 1
                            },
                            marker={
                                'size': 10,
                                'color': '#2e82bf'
                            },
                            showlegend=False
                        ) for _y, _name in zip(_arrays, _names)
                    ],
                    layout={
                        'xaxis': {
                            'title': _type.capitalize(),
                            'zeroline': False,
                            'tickmode': 'array',
                            'tickvals': _names,
                            'type': 'category'
                        },
                        'yaxis': {
                            'title': _name.capitalize() + ' network correlation',
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
                    _filename='plot_box_' + _type + '_low_con_' + str(_low_connectivity) + '_' + _name
                )

        # bar plot
        for _name, _sums in zip(['same', 'different'], [_same_highest, _different_highest]):
            if _name in _plots:
                _fig = go.Figure(
                    data=go.Bar(
                        x=_names,
                        y=[_type_sums[TIME_LAG_INDEX] / sum(_type_sums) for _type_sums in _sums],
                        marker={
                            'color': '#2e82bf'
                        }
                    ),
                    layout={
                        'xaxis': {
                            'title': _type.capitalize(),
                            'zeroline': False,
                            'tickmode': 'array',
                            'tickvals': _names,
                            'type': 'category'
                        },
                        'yaxis': {
                            'title': 'Lag ' + str(TIME_LAG) + ' highest correlation fraction',
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
                    _filename='plot_bar_' + _type + '_low_con_' + str(_low_connectivity) + '_' + _name
                )


if __name__ == '__main__':
    main()
