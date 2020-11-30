import os
from itertools import product

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0

PAIR_DISTANCE_RANGES = [(4, 6), (6, 8), (8, 10)]


def main(_band=True, _high_temporal_resolution=False):
    _y_arrays = [[] for _i in PAIR_DISTANCE_RANGES]
    _names_array = []
    for _distances_index, _distances_range in enumerate(PAIR_DISTANCE_RANGES):
        print('Pair distance range:', str(_distances_range))
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
        _tuples = filtering.by_pair_distance_range(_tuples, _distances_range)
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
                                                                   _keys=['experiment', 'series_id', 'group',
                                                                          'cell_id'])
        _fiber_densities = compute.fiber_densities(_windows_to_compute)

        _experiments_fiber_densities = {
            _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
            for _key in _windows_dictionary
        }

        _tuples_by_experiment = organize.by_experiment(_tuples)

        _higher_same_counter = 0
        _valid_tuples = []
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
                if len(_same_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_same_experiment):
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
                            if len(_same_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_different_experiment):
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
                                _y_arrays[_distances_index].append(_point_distance)
                                _higher_same_counter += 1
                            else:
                                _y_arrays[_distances_index].append(-_point_distance)

                            if _same_tuple not in _valid_tuples:
                                _valid_tuples.append(_same_tuple)

        print('Total tuples:', len(_valid_tuples))
        print('Total points:', len(_y_arrays[_distances_index]))
        print('Wilcoxon around the zero:')
        print(wilcoxon(_y_arrays[_distances_index]))
        print('Higher same amount:', _higher_same_counter / len(_y_arrays[_distances_index]))
        _names_array.append(str(_distances_range[0]) + '-' + str(_distances_range[1]))

    # plot
    _colors_array = config.colors(3)
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
                    'color': _color
                },
                showlegend=False
            ) for _y, _name, _color in zip(_y_arrays, _names_array, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Pair distance (cell diameter)',
                'zeroline': False,
                'type': 'category'
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
        _filename='plot_band_' + str(_band) + '_high_temporal_res_' + str(_high_temporal_resolution)
    )


if __name__ == '__main__':
    main()
