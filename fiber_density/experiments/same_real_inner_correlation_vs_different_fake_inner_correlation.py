import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
# based on time resolution
from plotting import save

EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
PAIR_DISTANCE_RANGE = [4, 10]
MINIMUM_CORRELATION_TIME_FRAMES = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}


def compute_fiber_densities(_offset_y=0.5, _high_time_resolution=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_real_fake_pairs(_experiments, _real_fake_pairs=False)
    _experiments_matched = organize.by_matched_real_and_fake(_experiments)
    print('Total matched pairs:', len(_experiments_matched))

    _arguments = []
    for _matched_tuple in _experiments_matched:
        for _tuple in _matched_tuple:
            _experiment, _series_id, _group = _tuple

            # stop when windows are overlapping
            _properties = load.group_properties(_experiment, _series_id, _group)
            _latest_time_frame = len(_properties['time_points'])
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
                    'offset_y': _offset_y,
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

    _tuples_by_experiment = organize.by_experiment(_experiments)

    # same (real, fake), different (real, fake)
    _correlations = [[[], []], [[], []]]
    _valid_real_tuples = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]
        _experiments_matched = organize.by_matched_real_and_fake(_experiment_tuples)
        print('Matched pairs:', len(_experiments_matched))

        for _same_index in tqdm(range(len(_experiments_matched)), desc='Main loop'):
            for _group_type_index in [0, 1]:
                _same_tuple = _experiments_matched[_same_index][_group_type_index]
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

                _same_left_cell_fiber_densities_filtered, _same_right_cell_fiber_densities_filtered = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _same_left_cell_fiber_densities, _same_right_cell_fiber_densities
                    )

                # ignore small arrays
                if len(_same_left_cell_fiber_densities_filtered) < \
                        MINIMUM_CORRELATION_TIME_FRAMES[_same_experiment]:
                    for _different_index in range(len(_experiments_matched)):
                        if _same_index != _different_index:
                            # for all combinations
                            for _i in range(4):
                                _correlations[0][_group_type_index].append(None)
                                _correlations[1][_group_type_index].append(None)
                    continue

                _same_correlation = compute_lib.correlation(
                    compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                    compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
                )
                for _different_index in range(len(_experiments_matched)):
                    if _same_index != _different_index:
                        _different_tuple = _experiments_matched[_different_index][_group_type_index]
                        _different_experiment, _different_series, _different_group = _different_tuple
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

                            _same_fiber_densities_filtered, _different_fiber_densities_filtered = \
                                compute.longest_same_indices_shared_in_borders_sub_array(
                                    _same_fiber_densities, _different_fiber_densities
                                )

                            # ignore small arrays
                            if len(_same_fiber_densities_filtered) < \
                                    MINIMUM_CORRELATION_TIME_FRAMES[_different_experiment]:
                                _correlations[0][_group_type_index].append(None)
                                _correlations[1][_group_type_index].append(None)
                                continue

                            _different_correlation = compute_lib.correlation(
                                compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                                compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                            )

                            _correlations[0][_group_type_index].append(_same_correlation)
                            _correlations[1][_group_type_index].append(_different_correlation)

                            if _group_type_index == 0 and _same_tuple not in _valid_real_tuples:
                                _valid_real_tuples.append(_same_tuple)

    print('Total real tuples:', len(_valid_real_tuples))
    _distances_from_y_equal_x = [[], []]
    _same_correlations, _different_correlations = _correlations
    _same_real_correlations, _same_fake_correlations = _same_correlations
    _different_real_correlations, _different_fake_correlations = _different_correlations
    for _same_real, _same_fake, _different_real, _different_fake in \
            zip(_same_real_correlations, _same_fake_correlations,
                _different_real_correlations, _different_fake_correlations):

        # one of the correlations is none - not valid
        if None in [_same_real, _same_fake, _different_real, _different_fake]:
            continue

        for _group_type_index, _same, _different in \
                zip([0, 1], [_same_real, _same_fake], [_different_real, _different_fake]):

            _point_distance = compute_lib.distance_from_a_point_to_a_line(
                _line=[-1, -1, 1, 1],
                _point=[_same, _different]
            )
            if _same > _different:
                _distances_from_y_equal_x[_group_type_index].append(_point_distance)
            else:
                _distances_from_y_equal_x[_group_type_index].append(-_point_distance)

    return _distances_from_y_equal_x


def main(_offset_y=0.5, _high_time_resolution=False):
    _distances_from_y_equal_x = compute_fiber_densities(_offset_y, _high_time_resolution)

    print('Total points:', len(_distances_from_y_equal_x[0]))
    print('Higher real same amount:', (np.array(_distances_from_y_equal_x[0]) > 0).sum() /
          len(_distances_from_y_equal_x[0]))
    print('Wilcoxon of real points:', wilcoxon(_distances_from_y_equal_x[0]))
    print('Higher fake same amount:', (np.array(_distances_from_y_equal_x[1]) > 0).sum() /
          len(_distances_from_y_equal_x[1]))
    print('Wilcoxon of fake points:', wilcoxon(_distances_from_y_equal_x[1]))
    _real_minus_fake = np.array(_distances_from_y_equal_x[0]) - np.array(_distances_from_y_equal_x[1])
    print('Real > fake amount:', (_real_minus_fake > 0).sum() /
          len(_real_minus_fake))
    print('Wilcoxon real & fake:', wilcoxon(_distances_from_y_equal_x[0], _distances_from_y_equal_x[1]))

    # box plot
    _colors_array = ['#844b00', '#ea8500']
    _names_array = ['Real', 'Fake']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y_array,
                name=_name,
                boxpoints=False,
                line={
                    'width': 1
                },
                marker={
                    'color': _color
                },
                showlegend=False
            ) for _y_array, _name, _color in zip(_distances_from_y_equal_x, _names_array, _colors_array)
        ],
        layout={
            'xaxis': {
                'zeroline': False
            },
            'yaxis': {
                'title': 'Same minus different correlation',
                'zeroline': False,
                'range': [-1, 1.1],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_box_high_time_res_' + str(_high_time_resolution) + '_offset_y_' + str(_offset_y)
    )

    # scatter plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_distances_from_y_equal_x[0],
            y=_distances_from_y_equal_x[1],
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Real pair',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Fake pair',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': -1,
                    'y1': 1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': 1,
                    'y1': -1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': 1,
                    'y1': 1,
                    'line': {
                        'color': 'red',
                        'width': 2
                    }
                }
            ]
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_res_' + str(_high_time_resolution) + '_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
