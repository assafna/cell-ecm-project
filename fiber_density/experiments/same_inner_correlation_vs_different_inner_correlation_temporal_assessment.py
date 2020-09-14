import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
TIME_FRAMES_STEPS = {
    'SN16': range(1, 6),
    'SN41': range(1, 18),
    'SN44': range(1, 18),
    'SN45': range(1, 18)
}
TIME_RESOLUTION = {
    'SN16': 15,
    'SN41': 5,
    'SN44': 5,
    'SN45': 5
}
OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
MINIMUM_CORRELATION_TIME_FRAMES = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}
GENERAL_MINIMUM_CORRELATION_TIME_FRAMES = {
    'SN16': 5,
    'SN41': 15,
    'SN44': 15,
    'SN45': 15
}


def main(_high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_real_pairs(_experiments, _real_pairs=REAL_CELLS)
    _experiments = filtering.by_fake_static_pairs(_experiments, _fake_static_pairs=STATIC)
    _experiments = filtering.by_band(_experiments)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
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

    _tuples_by_experiment = organize.by_experiment(_experiments)

    _y_arrays = [[] for _i in TIME_FRAMES_STEPS[EXPERIMENTS[_high_time_resolution][0]]]
    _x_array = []
    for _time_frame_index, _time_frame_every in enumerate(TIME_FRAMES_STEPS[EXPERIMENTS[_high_time_resolution][0]]):
        print('Time resolution (minutes):', _time_frame_every * TIME_RESOLUTION[EXPERIMENTS[_high_time_resolution][0]])
        _higher_same_counter = 0
        _valid_tuples = []
        for _time_frame_begin in range(_time_frame_every):
            for _experiment in _tuples_by_experiment:
                _experiment_tuples = _tuples_by_experiment[_experiment]

                for _same_index in range(len(_experiment_tuples)):
                    _same_tuple = _experiment_tuples[_same_index]
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
                        continue

                    _same_left_cell_fiber_densities_filtered = \
                        _same_left_cell_fiber_densities_filtered[_time_frame_begin::_time_frame_every]
                    _same_right_cell_fiber_densities_filtered = \
                        _same_right_cell_fiber_densities_filtered[_time_frame_begin::_time_frame_every]

                    # secondary ignore
                    if len(_same_left_cell_fiber_densities_filtered) < \
                            GENERAL_MINIMUM_CORRELATION_TIME_FRAMES[_experiment]:
                        continue

                    _same_correlation = compute_lib.correlation(
                        compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
                    )
                    for _different_index in range(len(_experiment_tuples)):
                        if _same_index != _different_index:
                            _different_tuple = _experiment_tuples[_different_index]
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

                                _same_fiber_densities_filtered, _different_fiber_densities_filtered = \
                                    compute.longest_same_indices_shared_in_borders_sub_array(
                                        _same_fiber_densities, _different_fiber_densities
                                    )

                                # ignore small arrays
                                if len(_same_fiber_densities_filtered) < \
                                        MINIMUM_CORRELATION_TIME_FRAMES[_different_experiment]:
                                    continue

                                _same_fiber_densities_filtered = \
                                    _same_fiber_densities_filtered[_time_frame_begin::_time_frame_every]
                                _different_fiber_densities_filtered = \
                                    _different_fiber_densities_filtered[_time_frame_begin::_time_frame_every]

                                # secondary ignore
                                if len(_same_fiber_densities_filtered) < \
                                        GENERAL_MINIMUM_CORRELATION_TIME_FRAMES[_experiment]:
                                    continue

                                _different_correlation = compute_lib.correlation(
                                    compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                                    compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                                )

                                _point_distance = compute_lib.distance_from_a_point_to_a_line(
                                    _line=[-1, -1, 1, 1],
                                    _point=[_same_correlation, _different_correlation]
                                )
                                if _same_correlation > _different_correlation:
                                    _y_arrays[_time_frame_index].append(_point_distance)
                                    _higher_same_counter += 1
                                else:
                                    _y_arrays[_time_frame_index].append(-_point_distance)

                                if _same_tuple not in _valid_tuples:
                                    _valid_tuples.append(_same_tuple)

        print('Total tuples:', len(_valid_tuples))
        print('Total points:', len(_y_arrays[_time_frame_index]))
        print('Wilcoxon around the zero:')
        print(wilcoxon(_y_arrays[_time_frame_index]))
        print('Higher same amount:', _higher_same_counter / len(_y_arrays[_time_frame_index]))
        _x_array.append(_time_frame_every * TIME_RESOLUTION[EXPERIMENTS[_high_time_resolution][0]])

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_x_array,
            y=[np.mean(_array) for _array in _y_arrays],
            error_y={
                'type': 'data',
                'array': [np.std(_array) for _array in _y_arrays],
                'thickness': 1,
                'color': '#ea8500'
            },
            mode='markers',
            marker={
                'size': 15,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Time resolution (minutes)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Same minus different correlation',
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
        _filename='plot_high_time_res_' + str(_high_time_resolution)
    )


if __name__ == '__main__':
    main()
