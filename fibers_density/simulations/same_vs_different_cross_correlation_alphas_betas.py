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


def main(_type='alpha', _low_connectivity=False, _plots=None, _plot_types=None):
    if _plots is None:
        _plots = ['same', 'different']
    if _plot_types is None:
        _plot_types = ['stacked_bar', 'box', 'bar']

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

        # stacked bar plot
        if 'stacked_bar' in _plot_types:
            for _name, _sums in zip(['same', 'different'], [_same_highest, _different_highest]):
                if _name in _plots:

                    _y_arrays = [[], [], []]
                    for _type_sums in _sums:
                        _left_wins, _none_wins, _right_wins = 0, 0, 0
                        for _time_lag, _type_sum in zip(same_vs_different_cross_correlation.TIME_LAGS, _type_sums):
                            if _time_lag > 0:
                                _left_wins += _type_sum
                            elif _time_lag < 0:
                                _right_wins += _type_sum
                            else:
                                _none_wins += _type_sum
                        _total = sum(_type_sums)
                        _y_arrays[0].append(_left_wins / _total)
                        _y_arrays[1].append(_none_wins / _total)
                        _y_arrays[2].append(_right_wins / _total)

                    _colors_array = ['#011f4b', '#005b96', '#74c2e8']
                    _fig = go.Figure(
                        data=[
                            go.Bar(
                                x=_names,
                                y=_y_array,
                                name=_name,
                                marker={
                                    'color': _color
                                }
                            ) for _name, _y_array, _color in
                            zip(['Left cell', 'None', 'Right cell'], _y_arrays, _colors_array)
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
                                'title': 'Highest correlation fraction',
                                'range': [0, 1.1],
                                'zeroline': False,
                                'tickmode': 'array',
                                'tickvals': [0, 0.5, 1]
                            },
                            'barmode': 'stack',
                            'legend': {
                                'bordercolor': 'black',
                                'borderwidth': 2,
                                'bgcolor': 'white'
                            },
                        }
                    )

                    save.to_html(
                        _fig=_fig,
                        _path=os.path.join(paths.PLOTS, save.get_module_name()),
                        _filename='plot_stacked_bar_' + _type + '_low_con_' + str(_low_connectivity) + '_' + _name
                    )

        # box plot
        if 'box' in _plot_types:
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
                                'title': 'Inner correlation' if _name == 'same' else 'Different network correlation',
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
        if 'bar' in _plot_types:
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
