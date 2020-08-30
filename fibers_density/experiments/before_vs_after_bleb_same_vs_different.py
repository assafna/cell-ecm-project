import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = ['SN26_BlebAdded']
AFTER_BLEB_INJECTION_FIRST_TIME_POINT = 9
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
DIRECTION = 'inside'
MINIMUM_TIME_POINTS = 20
MINIMUM_TIME_POINTS_CORRELATION = 8
DERIVATIVE = 1


def main(_band=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_time_points_amount(_experiments, _time_points=MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance_range(_experiments, _distance_range=CELLS_DISTANCE_RANGE)
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

    _tuples_by_experiment = organize.by_experiment(_experiments)

    # same (before, after), different (before, after)
    _correlations = [[[], []], [[], []]]
    _valid_real_tuples = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _same_index in tqdm(range(len(_experiment_tuples)), desc='Main loop'):
            _same_tuple = _experiment_tuples[_same_index]
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

            _same_before_left_cell_fiber_densities = \
                _same_left_cell_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]
            _same_before_right_cell_fiber_densities = \
                _same_right_cell_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]

            _same_after_left_cell_fiber_densities = \
                _same_left_cell_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]
            _same_after_right_cell_fiber_densities = \
                _same_right_cell_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]

            _same_before_left_cell_fibers_densities_filtered, _same_before_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_before_left_cell_fiber_densities, _same_before_right_cell_fiber_densities
                )

            _same_after_left_cell_fibers_densities_filtered, _same_after_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_after_left_cell_fiber_densities, _same_after_right_cell_fiber_densities
                )

            # ignore small arrays
            if len(_same_before_left_cell_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION or \
                    len(_same_after_left_cell_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION:
                continue

            _same_before_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_before_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_before_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            _same_after_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_after_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_after_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )

            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
                    _different_experiment, _different_series, _different_group = _different_tuple
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

                        _same_before_fiber_densities = \
                            _same_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]
                        _same_after_fiber_densities = \
                            _same_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]

                        _different_before_fiber_densities = \
                            _different_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]
                        _different_after_fiber_densities = \
                            _different_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]

                        _same_before_fibers_densities_filtered, _different_before_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_before_fiber_densities, _different_before_fiber_densities
                            )

                        _same_after_fibers_densities_filtered, _different_after_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_after_fiber_densities, _different_after_fiber_densities
                            )

                        # ignore small arrays
                        if len(_same_before_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION or \
                                len(_same_after_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION:
                            continue

                        _different_before_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_before_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_before_fibers_densities_filtered, _n=DERIVATIVE)
                        )
                        _different_after_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_after_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_after_fibers_densities_filtered, _n=DERIVATIVE)
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
    _colors_array = ['#844b00', '#ea8500']
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
