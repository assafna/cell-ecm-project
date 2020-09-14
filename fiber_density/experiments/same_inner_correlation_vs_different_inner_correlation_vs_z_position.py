import os
from itertools import product

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

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
OFFSET_X = 0
OFFSET_Y = 0
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


def compute_fiber_densities(_band=True, _high_time_resolution=False):
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

    _distances_from_y_equal_x = []
    _z_positions_array = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _same_index in tqdm(range(len(_experiment_tuples)), desc='Main loop'):
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

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )

            _same_group_mean_z_position = \
                compute.group_mean_z_position_from_substrate(_same_experiment, _same_series, _same_group)

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

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                        )

                        _point_distance = compute_lib.distance_from_a_point_to_a_line(
                            _line=[-1, -1, 1, 1],
                            _point=[_same_correlation, _different_correlation]
                        )
                        if _same_correlation > _different_correlation:
                            _distances_from_y_equal_x.append(_point_distance)
                        else:
                            _distances_from_y_equal_x.append(-_point_distance)

                        _z_positions_array.append(_same_group_mean_z_position)

    print('Total points:', len(_distances_from_y_equal_x))
    print('Wilcoxon of distances from y = x around the zero:')
    print(wilcoxon(_distances_from_y_equal_x))
    print('Pearson correlation of distances from y = x and z position distances:')
    print(compute_lib.correlation(_distances_from_y_equal_x, _z_positions_array, _with_p_value=True))

    return _distances_from_y_equal_x, _z_positions_array


def main(_band=True, _high_time_resolution=False):
    _distances_from_y_equal_x, _z_positions_array = compute_fiber_densities(_band, _high_time_resolution)

    _min_z_position, _max_z_position = min(_z_positions_array), max(_z_positions_array)
    _z_positions_array_normalized = [(_z_position_value - _min_z_position) / (_max_z_position - _min_z_position)
                                     for _z_position_value in _z_positions_array]

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_z_positions_array_normalized,
            y=_distances_from_y_equal_x,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Normalized mean Z distance',
                'zeroline': False,
                'range': [-0.1, 1.2],
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Distance from y = x',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_band_' + str(_band) + '_high_time_res_' + str(_high_time_resolution)
    )


if __name__ == '__main__':
    main()
