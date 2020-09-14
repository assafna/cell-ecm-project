import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
TIME_RESOLUTION = {
    False: 15,
    True: 5
}
OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
DIRECTION = 'inside'
MINIMUM_CORRELATION_TIME_FRAMES = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}
TIME_LAGS = {
    False: [0, 1, 2],
    True: [0, 1, 2, 3, 4, 5]
}


def compute_fiber_densities(_band=True, _high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_real_pairs(_experiments, _real_pairs=REAL_CELLS)
    _experiments = filtering.by_fake_static_pairs(_experiments, _fake_static_pairs=STATIC)
    _experiments = filtering.by_band(_experiments, _band=_band)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_frame = len(_properties['time_points'])
        if DIRECTION == 'inside':
            for _time_frame in range(len(_properties['time_points'])):
                _pair_distance = \
                    compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame)
                if _pair_distance - 1 - OFFSET_X * 2 < QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER * 2:
                    _latest_time_frame = _time_frame - 1
                    break

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _same_correlation_vs_time_lag = {}
    _same_time_lags_arrays = [[] for _i in TIME_LAGS[_high_time_resolution]]
    _different_time_lags_arrays = [[] for _i in TIME_LAGS[_high_time_resolution]]
    _same_time_lags_highest = [0 for _i in TIME_LAGS[_high_time_resolution]]
    _different_time_lags_highest = [0 for _i in TIME_LAGS[_high_time_resolution]]
    _valid_tuples = []
    for _same_index in tqdm(range(len(_experiments)), desc='Main loop'):
        _same_tuple = _experiments[_same_index]
        _same_experiment, _same_series, _same_group = _same_tuple

        _same_left_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_same_experiment, _same_series, _same_group, 'left_cell')
            ]
        _same_right_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_same_experiment, _same_series, _same_group, 'right_cell')
            ]

        _same_properties = \
            load.group_properties(_same_experiment, _same_series, _same_group)
        _same_left_cell_fiber_densities = compute.remove_blacklist(
            _same_experiment,
            _same_series,
            _same_properties['cells_ids']['left_cell'],
            _same_left_cell_fiber_densities
        )
        _same_right_cell_fiber_densities = compute.remove_blacklist(
            _same_experiment,
            _same_series,
            _same_properties['cells_ids']['right_cell'],
            _same_right_cell_fiber_densities
        )

        # time lag
        _same_highest_correlation = -1.1
        _same_highest_correlation_time_lag_index = 0
        _same_correlation_vs_time_lag[_same_tuple] = []
        for _time_lag_index, _time_lag in enumerate(TIME_LAGS[_high_time_resolution]):

            # choose either negative or positive lag
            for _symbol in [-1, 1]:

                # if no time lag consider it only once
                if _time_lag == 0 and _symbol == -1:
                    continue

                _time_lag_symbol = _time_lag * _symbol

                if _time_lag_symbol > 0:
                    _same_left_cell_fiber_densities_time_lag = _same_left_cell_fiber_densities[:-_time_lag_symbol]
                    _same_right_cell_fiber_densities_time_lag = _same_right_cell_fiber_densities[_time_lag_symbol:]
                elif _time_lag_symbol < 0:
                    _same_left_cell_fiber_densities_time_lag = _same_left_cell_fiber_densities[-_time_lag_symbol:]
                    _same_right_cell_fiber_densities_time_lag = _same_right_cell_fiber_densities[:_time_lag_symbol]
                else:
                    _same_left_cell_fiber_densities_time_lag = _same_left_cell_fiber_densities
                    _same_right_cell_fiber_densities_time_lag = _same_right_cell_fiber_densities

                _same_left_cell_fiber_densities_filtered, _same_right_cell_fiber_densities_filtered = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _same_left_cell_fiber_densities_time_lag, _same_right_cell_fiber_densities_time_lag
                    )

                # ignore small arrays
                if len(_same_left_cell_fiber_densities_filtered) < \
                        MINIMUM_CORRELATION_TIME_FRAMES[_same_experiment]:
                    _same_correlation_vs_time_lag[_same_tuple].append(None)
                    continue

                _same_correlation = compute_lib.correlation(
                    compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                    compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
                )

                _same_time_lags_arrays[_time_lag_index].append(_same_correlation)
                _same_correlation_vs_time_lag[_same_tuple].append(_same_correlation)

                if _same_correlation > _same_highest_correlation:
                    _same_highest_correlation = _same_correlation
                    _same_highest_correlation_time_lag_index = _time_lag_index

        _same_time_lags_highest[_same_highest_correlation_time_lag_index] += 1

        for _different_index in range(len(_experiments)):
            if _same_index != _different_index:
                _different_tuple = _experiments[_different_index]
                _different_experiment, _different_series, _different_group = \
                    _different_tuple
                for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                 ['left_cell', 'right_cell']):
                    _same_fiber_densities = _experiments_fiber_densities[(
                        _same_experiment,
                        _same_series,
                        _same_group,
                        _same_cell_id
                    )]
                    _different_fiber_densities = _experiments_fiber_densities[(
                        _different_experiment,
                        _different_series,
                        _different_group,
                        _different_cell_id
                    )]

                    _different_properties = load.group_properties(
                        _different_experiment, _different_series, _different_group
                    )
                    _same_fiber_densities = compute.remove_blacklist(
                        _same_experiment,
                        _same_series,
                        _same_properties['cells_ids'][_same_cell_id],
                        _same_fiber_densities
                    )
                    _different_fiber_densities = compute.remove_blacklist(
                        _different_experiment,
                        _different_series,
                        _different_properties['cells_ids'][_different_cell_id],
                        _different_fiber_densities
                    )

                    # time lag
                    _different_highest_correlation = -1.1
                    _different_highest_correlation_time_lag_index = 0
                    for _time_lag_index, _time_lag in enumerate(TIME_LAGS[_high_time_resolution]):

                        # choose either negative or positive lag
                        for _symbol in [-1, 1]:

                            # if no time lag consider it only once
                            if _time_lag == 0 and _symbol == -1:
                                continue

                            _time_lag_symbol = _time_lag * _symbol

                            if _time_lag_symbol > 0:
                                _same_fiber_densities_time_lag = _same_fiber_densities[:-_time_lag_symbol]
                                _different_fiber_densities_time_lag = _different_fiber_densities[_time_lag_symbol:]
                            elif _time_lag_symbol < 0:
                                _same_fiber_densities_time_lag = _same_fiber_densities[-_time_lag_symbol:]
                                _different_fiber_densities_time_lag = _different_fiber_densities[:_time_lag_symbol]
                            else:
                                _same_fiber_densities_time_lag = _same_fiber_densities
                                _different_fiber_densities_time_lag = _different_fiber_densities

                            _same_fiber_densities_filtered, _different_fiber_densities_filtered = \
                                compute.longest_same_indices_shared_in_borders_sub_array(
                                    _same_fiber_densities_time_lag, _different_fiber_densities_time_lag
                                )

                            # ignore small arrays
                            if len(_same_fiber_densities_filtered) < \
                                    MINIMUM_CORRELATION_TIME_FRAMES[_different_experiment]:
                                continue

                            _different_correlation = compute_lib.correlation(
                                compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                                compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                            )

                            _different_time_lags_arrays[_time_lag_index].append(_different_correlation)

                            if _different_correlation > _different_highest_correlation:
                                _different_highest_correlation = _different_correlation
                                _different_highest_correlation_time_lag_index = _time_lag_index

                            if _same_tuple not in _valid_tuples:
                                _valid_tuples.append(_same_tuple)

                    _different_time_lags_highest[_different_highest_correlation_time_lag_index] += 1

    print('Total tuples:', len(_valid_tuples))

    return _same_correlation_vs_time_lag, _same_time_lags_arrays, _different_time_lags_arrays, \
        _same_time_lags_highest, _different_time_lags_highest


def main(_band=True, _high_time_resolution=True, _plots=None, _plot_types=None):
    if _plots is None:
        _plots = ['same', 'different']
    if _plot_types is None:
        _plot_types = ['scatter', 'box', 'bar']

    _same_correlation_vs_time_lag, _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, \
        _different_time_lags_highest = compute_fiber_densities(_band, _high_time_resolution)

    if _plots is not None:

        # individual plots
        if 'scatter' in _plot_types:
            for _same_tuple in _same_correlation_vs_time_lag:
                _experiment, _series_id, _group = _same_tuple
                _fig = go.Figure(
                    data=go.Scatter(
                        x=np.array(TIME_LAGS[_high_time_resolution]) * TIME_RESOLUTION[_high_time_resolution],
                        y=_same_correlation_vs_time_lag[_same_tuple],
                        mode='markers',
                        marker={
                            'size': 25,
                            'color': '#ea8500'
                        }
                    ),
                    layout={
                        'xaxis': {
                            'title': 'Time lag (minutes)',
                            'zeroline': False,
                            'tickmode': 'array',
                            'tickvals': np.array(TIME_LAGS[_high_time_resolution]) *
                                        TIME_RESOLUTION[_high_time_resolution]
                        },
                        'yaxis': {
                            'title': 'Correlation',
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
                    _filename='plot_' + _experiment + '_' + str(_series_id) + '_' + _group
                )

        # box plots
        if 'box' in _plot_types:
            for _name, _arrays in zip(['same', 'different'], [_same_time_lags_arrays, _different_time_lags_arrays]):
                if _name in _plots:
                    _fig = go.Figure(
                        data=[
                            go.Box(
                                y=_y,
                                name=_time_lag * TIME_RESOLUTION[_high_time_resolution],
                                boxpoints=False,
                                line={
                                    'width': 1
                                },
                                marker={
                                    'size': 10,
                                    'color': '#ea8500'
                                },
                                showlegend=False
                            ) for _y, _time_lag in zip(_arrays, TIME_LAGS[_high_time_resolution])
                        ],
                        layout={
                            'xaxis': {
                                'title': 'Time lag (minutes)',
                                'zeroline': False,
                                'tickmode': 'array',
                                'tickvals': np.array(TIME_LAGS[_high_time_resolution]) *
                                            TIME_RESOLUTION[_high_time_resolution]
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
                        _filename='plot_box_high_time_res_' + str(_high_time_resolution) + '_' + _name
                    )

        # bar plot
        if 'bar' in _plot_types:
            for _name, _sums in zip(['same', 'different'], [_same_time_lags_highest, _different_time_lags_highest]):
                if _name in _plots:
                    _fig = go.Figure(
                        data=go.Bar(
                            x=np.array(TIME_LAGS[_high_time_resolution]) * TIME_RESOLUTION[_high_time_resolution],
                            y=np.array(_sums) / sum(_sums),
                            marker={
                                'color': '#ea8500'
                            }
                        ),
                        layout={
                            'xaxis': {
                                'title': 'Time lag (minutes)',
                                'zeroline': False,
                                'tickmode': 'array',
                                'tickvals': np.array(TIME_LAGS[_high_time_resolution]) *
                                            TIME_RESOLUTION[_high_time_resolution]
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
                        _filename='plot_bar_high_time_res_' + str(_high_time_resolution) + '_' + _name
                    )


if __name__ == '__main__':
    main()
