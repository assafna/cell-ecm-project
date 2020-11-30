import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]

TIME_LAGS = {
    False: [0, 1, 2],
    True: [0, 1, 2, 3, 4, 5]
}


def compute_fiber_densities(_band=True, _high_temporal_resolution=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
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
                'direction': 'inside',
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
    _same_time_lags_arrays = [[] for _i in TIME_LAGS[_high_temporal_resolution]]
    _different_time_lags_arrays = [[] for _i in TIME_LAGS[_high_temporal_resolution]]
    _same_time_lags_highest = [0 for _i in TIME_LAGS[_high_temporal_resolution]]
    _different_time_lags_highest = [0 for _i in TIME_LAGS[_high_temporal_resolution]]
    _valid_tuples = []
    for _same_index in tqdm(range(len(_tuples)), desc='Main loop'):
        _same_tuple = _tuples[_same_index]
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
        for _time_lag_index, _time_lag in enumerate(TIME_LAGS[_high_temporal_resolution]):

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
                if len(_same_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_same_experiment):
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

        for _different_index in range(len(_tuples)):
            if _same_index != _different_index:
                _different_tuple = _tuples[_different_index]
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
                    for _time_lag_index, _time_lag in enumerate(TIME_LAGS[_high_temporal_resolution]):

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
                            if len(_same_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_different_experiment):
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


def main(_band=True, _high_temporal_resolution=True, _plots=None, _plot_types=None):
    if _plots is None:
        _plots = ['same', 'different']
    if _plot_types is None:
        _plot_types = ['scatter', 'box', 'bar']

    _same_correlation_vs_time_lag, _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, \
        _different_time_lags_highest = compute_fiber_densities(_band, _high_temporal_resolution)

    if _plots is not None:

        # individual plots
        if 'scatter' in _plot_types:
            for _same_tuple in _same_correlation_vs_time_lag:
                _experiment, _series_id, _group = _same_tuple
                _temporal_resolution = compute.temporal_resolution_in_minutes(_experiment)
                _fig = go.Figure(
                    data=go.Scatter(
                        x=np.array(TIME_LAGS[_high_temporal_resolution]) * _temporal_resolution,
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
                            'tickvals': np.array(TIME_LAGS[_high_temporal_resolution]) * _temporal_resolution
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
                                name=_time_lag * 5 if _high_temporal_resolution else 15,
                                boxpoints=False,
                                line={
                                    'width': 1
                                },
                                marker={
                                    'size': 10,
                                    'color': '#ea8500'
                                },
                                showlegend=False
                            ) for _y, _time_lag in zip(_arrays, TIME_LAGS[_high_temporal_resolution])
                        ],
                        layout={
                            'xaxis': {
                                'title': 'Time lag (minutes)',
                                'zeroline': False,
                                'tickmode': 'array',
                                'tickvals': np.array(TIME_LAGS[_high_temporal_resolution]) * 5 if _high_temporal_resolution else 15
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
                        _filename='plot_box_high_temporal_res_' + str(_high_temporal_resolution) + '_' + _name
                    )

        # bar plot
        if 'bar' in _plot_types:
            for _name, _sums in zip(['same', 'different'], [_same_time_lags_highest, _different_time_lags_highest]):
                if _name in _plots:
                    _fig = go.Figure(
                        data=go.Bar(
                            x=np.array(TIME_LAGS[_high_temporal_resolution]) * 5 if _high_temporal_resolution else 15,
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
                                'tickvals': np.array(TIME_LAGS[_high_temporal_resolution]) * 5 if _high_temporal_resolution else 15
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
                        _filename='plot_bar_high_temporal_res_' + str(_high_temporal_resolution) + '_' + _name
                    )


if __name__ == '__main__':
    main()
