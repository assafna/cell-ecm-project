import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41']
}
TIME_RESOLUTION = {
    False: 15,
    True: 5
}
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
DIRECTION = 'inside'
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 10,
    'SN18': 10,
    'SN41': 40,
    'SN44': 40
}
TIME_LAGS = [-2, -1, 0, 1, 2]


def compute_fibers_densities(_band=True, _high_time_resolution=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=_band)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
        if DIRECTION == 'inside':
            for _time_point in range(len(_properties['time_points'])):
                _cells_distance = \
                    compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
                if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                    _latest_time_point = _time_point - 1
                    break

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _same_time_lags_arrays = [[] for _i in TIME_LAGS]
    _different_time_lags_arrays = [[] for _i in TIME_LAGS]
    _same_time_lags_highest = [0 for _i in TIME_LAGS]
    _different_time_lags_highest = [0 for _i in TIME_LAGS]
    for _same_index in tqdm(range(len(_experiments)), desc='Main loop'):
        _same_tuple = _experiments[_same_index]
        _same_experiment, _same_series, _same_group = _same_tuple

        _same_left_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_same_experiment, _same_series, _same_group, 'left_cell')
            ]
        _same_right_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_same_experiment, _same_series, _same_group, 'right_cell')
            ]

        _same_properties = \
            load.group_properties(_same_experiment, _same_series, _same_group)
        _same_left_cell_fibers_densities = compute.remove_blacklist(
            _same_experiment,
            _same_series,
            _same_properties['cells_ids']['left_cell'],
            _same_left_cell_fibers_densities
        )
        _same_right_cell_fibers_densities = compute.remove_blacklist(
            _same_experiment,
            _same_series,
            _same_properties['cells_ids']['right_cell'],
            _same_right_cell_fibers_densities
        )

        # time lag
        _same_highest_correlation = -1.1
        _same_highest_correlation_time_lag_index = 0
        for _time_lag_index, _time_lag in enumerate(TIME_LAGS):
            if _time_lag > 0:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities[:-_time_lag]
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities[_time_lag:]
            elif _time_lag < 0:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities[-_time_lag:]
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities[:_time_lag]
            else:
                _same_left_cell_fibers_densities_time_lag = _same_left_cell_fibers_densities
                _same_right_cell_fibers_densities_time_lag = _same_right_cell_fibers_densities

            _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fibers_densities_time_lag, _same_right_cell_fibers_densities_time_lag
                )

            # ignore small arrays
            if len(_same_left_cell_fibers_densities_filtered) < \
                    MINIMUM_CORRELATION_TIME_POINTS[_same_experiment]:
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )

            _same_time_lags_arrays[_time_lag_index].append(_same_correlation)

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
                    _same_fibers_densities = _experiments_fibers_densities[(
                        _same_experiment,
                        _same_series,
                        _same_group,
                        _same_cell_id
                    )]
                    _different_fibers_densities = _experiments_fibers_densities[(
                        _different_experiment,
                        _different_series,
                        _different_group,
                        _different_cell_id
                    )]

                    _different_properties = load.group_properties(
                        _different_experiment, _different_series, _different_group
                    )
                    _same_fibers_densities = compute.remove_blacklist(
                        _same_experiment,
                        _same_series,
                        _same_properties['cells_ids'][_same_cell_id],
                        _same_fibers_densities
                    )
                    _different_fibers_densities = compute.remove_blacklist(
                        _different_experiment,
                        _different_series,
                        _different_properties['cells_ids'][_different_cell_id],
                        _different_fibers_densities
                    )

                    # time lag
                    _different_highest_correlation = -1.1
                    _different_highest_correlation_time_lag_index = 0
                    for _time_lag_index, _time_lag in enumerate(TIME_LAGS):
                        if _time_lag > 0:
                            _same_fibers_densities_time_lag = _same_fibers_densities[:-_time_lag]
                            _different_fibers_densities_time_lag = _different_fibers_densities[_time_lag:]
                        elif _time_lag < 0:
                            _same_fibers_densities_time_lag = _same_fibers_densities[-_time_lag:]
                            _different_fibers_densities_time_lag = _different_fibers_densities[:_time_lag]
                        else:
                            _same_fibers_densities_time_lag = _same_fibers_densities
                            _different_fibers_densities_time_lag = _different_fibers_densities

                        _same_fibers_densities_filtered, _different_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fibers_densities_time_lag, _different_fibers_densities_time_lag
                            )

                        # ignore small arrays
                        if len(_same_fibers_densities_filtered) < \
                                MINIMUM_CORRELATION_TIME_POINTS[_different_experiment]:
                            continue

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
                        )

                        _different_time_lags_arrays[_time_lag_index].append(_different_correlation)

                        if _different_correlation > _different_highest_correlation:
                            _different_highest_correlation = _different_correlation
                            _different_highest_correlation_time_lag_index = _time_lag_index

                    _different_time_lags_highest[_different_highest_correlation_time_lag_index] += 1

    return _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, _different_time_lags_highest


def main(_band=True, _high_time_resolution=False):
    _same_time_lags_arrays, _different_time_lags_arrays, _same_time_lags_highest, _different_time_lags_highest = \
        compute_fibers_densities(_band, _high_time_resolution)

    _colors_array = ['#844b00', '#ea8500', '#edbc80', '#ea8500', '#844b00']

    # box plots
    for _name, _arrays in zip(['same', 'different'], [_same_time_lags_arrays, _different_time_lags_arrays]):
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
                        'color': _color
                    },
                    showlegend=False
                ) for _y, _time_lag, _color in zip(_arrays, TIME_LAGS, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Time lag (minutes)',
                    'zeroline': False
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
            _filename='plot_box_high_time_res_' + str(_high_time_resolution) + '_' + _name
        )

    # bar plot
    for _name, _sums in zip(['same', 'different'], [_same_time_lags_highest, _different_time_lags_highest]):
        _fig = go.Figure(
            data=go.Bar(
                x=np.array(TIME_LAGS) * TIME_RESOLUTION[_high_time_resolution],
                y=np.array(_sums) / sum(_sums),
                marker={
                    'color': _colors_array
                }
            ),
            layout={
                'xaxis': {
                    'title': 'Time lag (minutes)',
                    'zeroline': False
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
