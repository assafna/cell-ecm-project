import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    AFTER_BLEB_INJECTION_FIRST_TIME_FRAME, DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False


def main(_band=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=None,
        _is_bleb=True,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, _distance_range=PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=REAL_CELLS)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=STATIC)
    _tuples = filtering.by_band(_tuples, _band=_band)
    _tuples = filtering.by_bleb_from_start(_tuples, _from_start=False)
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
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _tuples_by_experiment = organize.by_experiment(_tuples)

    # same (before, after), different (before, after)
    _correlations = [[[], []], [[], []]]
    _valid_real_tuples = []
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

            _same_before_left_cell_fiber_densities = \
                _same_left_cell_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]]
            _same_before_right_cell_fiber_densities = \
                _same_right_cell_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]]

            _same_after_left_cell_fiber_densities = \
                _same_left_cell_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]:]
            _same_after_right_cell_fiber_densities = \
                _same_right_cell_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]:]

            _same_before_left_cell_fiber_densities_filtered, _same_before_right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_before_left_cell_fiber_densities, _same_before_right_cell_fiber_densities
                )

            _same_after_left_cell_fiber_densities_filtered, _same_after_right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_after_left_cell_fiber_densities, _same_after_right_cell_fiber_densities
                )

            # ignore small arrays
            _minimum_time_frame_for_correlation = compute.minimum_time_frames_for_correlation(_same_experiment)
            if len(_same_before_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation or \
                    len(_same_after_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation:
                continue

            _same_before_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_before_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_before_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )
            _same_after_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_after_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_after_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )

            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
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

                        _same_before_fiber_densities = \
                            _same_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]]
                        _same_after_fiber_densities = \
                            _same_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_same_experiment]:]

                        _different_before_fiber_densities = \
                            _different_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_different_experiment]]
                        _different_after_fiber_densities = \
                            _different_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_different_experiment]:]

                        _same_before_fiber_densities_filtered, _different_before_fiber_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_before_fiber_densities, _different_before_fiber_densities
                            )

                        _same_after_fiber_densities_filtered, _different_after_fiber_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_after_fiber_densities, _different_after_fiber_densities
                            )

                        # ignore small arrays
                        if len(_same_before_fiber_densities_filtered) < _minimum_time_frame_for_correlation or \
                                len(_same_after_fiber_densities_filtered) < _minimum_time_frame_for_correlation:
                            continue

                        _different_before_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_before_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_before_fiber_densities_filtered, _n=DERIVATIVE)
                        )
                        _different_after_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_after_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_after_fiber_densities_filtered, _n=DERIVATIVE)
                        )

                        _correlations[0][0].append(_same_before_correlation)
                        _correlations[0][1].append(_same_after_correlation)
                        _correlations[1][0].append(_different_before_correlation)
                        _correlations[1][1].append(_different_after_correlation)

                        if _same_tuple not in _valid_real_tuples:
                            _valid_real_tuples.append(_same_tuple)

    print('Total tuples:', len(_valid_real_tuples))
    _distances_from_y_equal_x = [[], []]
    _same_correlations, _different_correlations = _correlations
    _same_before_correlations, _same_after_correlations = _same_correlations
    _different_before_correlations, _different_after_correlations = _different_correlations
    for _same_before, _same_after, _different_before, _different_after in \
            zip(_same_before_correlations, _same_after_correlations,
                _different_before_correlations, _different_after_correlations):

        for _group_type_index, _same, _different in \
                zip([0, 1], [_same_before, _same_after], [_different_before, _different_after]):

            _point_distance = compute_lib.distance_from_a_point_to_a_line(
                _line=[-1, -1, 1, 1],
                _point=[_same, _different]
            )
            if _same > _different:
                _distances_from_y_equal_x[_group_type_index].append(_point_distance)
            else:
                _distances_from_y_equal_x[_group_type_index].append(-_point_distance)

    print('Total points:', len(_distances_from_y_equal_x[0]))
    print('Higher before same amount:', (np.array(_distances_from_y_equal_x[0]) > 0).sum() /
          len(_distances_from_y_equal_x[0]))
    print('Wilcoxon of before points:', wilcoxon(_distances_from_y_equal_x[0]))
    print('Higher after same amount:', (np.array(_distances_from_y_equal_x[1]) > 0).sum() /
          len(_distances_from_y_equal_x[1]))
    print('Wilcoxon of after points:', wilcoxon(_distances_from_y_equal_x[1]))
    _before_minus_after = np.array(_distances_from_y_equal_x[0]) - np.array(_distances_from_y_equal_x[1])
    print('Before > after amount:', (_before_minus_after > 0).sum() /
          len(_before_minus_after))
    print('Wilcoxon before & after:', wilcoxon(_distances_from_y_equal_x[0], _distances_from_y_equal_x[1]))

    # box plot
    _colors_array = config.colors(2)
    _names_array = ['Before', 'After']
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
        _filename='plot_box'
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
                'title': 'Before bleb',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'After bleb',
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
