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
EXPERIMENT = 'SN41'
TIME_POINTS_STEPS = {
    'SN16': range(1, 6),
    'SN41': range(1, 18)
}
TIME_RESOLUTION = {
    'SN16': 15,
    'SN41': 5
}
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}
GENERAL_MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 5,
    'SN41': 15
}


def main():
    _experiments = load.experiment_groups_as_tuples(EXPERIMENT)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
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
                'direction': 'inside',
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _y_arrays = [[] for _i in TIME_POINTS_STEPS[EXPERIMENT]]
    _x_array = []
    for _time_point_index, _time_point_every in enumerate(TIME_POINTS_STEPS[EXPERIMENT]):
        print('Time resolution (minutes):', _time_point_every * TIME_RESOLUTION[EXPERIMENT])
        _higher_same_counter = 0
        for _time_point_begin in range(_time_point_every):
            for _same_index in range(len(_experiments)):
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

                _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _same_left_cell_fibers_densities, _same_right_cell_fibers_densities
                    )

                # ignore small arrays
                if len(_same_left_cell_fibers_densities_filtered) < \
                        MINIMUM_CORRELATION_TIME_POINTS[_same_experiment]:
                    continue

                _same_left_cell_fibers_densities_filtered = \
                    _same_left_cell_fibers_densities_filtered[_time_point_begin::_time_point_every]
                _same_right_cell_fibers_densities_filtered = \
                    _same_right_cell_fibers_densities_filtered[_time_point_begin::_time_point_every]

                # secondary ignore
                if len(_same_left_cell_fibers_densities_filtered) < GENERAL_MINIMUM_CORRELATION_TIME_POINTS[EXPERIMENT]:
                    continue

                _same_correlation = compute_lib.correlation(
                    compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                    compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
                )
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

                            _same_fibers_densities_filtered, _different_fibers_densities_filtered = \
                                compute.longest_same_indices_shared_in_borders_sub_array(
                                    _same_fibers_densities, _different_fibers_densities
                                )

                            # ignore small arrays
                            if len(_same_fibers_densities_filtered) < \
                                    MINIMUM_CORRELATION_TIME_POINTS[_different_experiment]:
                                continue

                            _same_fibers_densities_filtered = \
                                _same_fibers_densities_filtered[_time_point_begin::_time_point_every]
                            _different_fibers_densities_filtered = \
                                _different_fibers_densities_filtered[_time_point_begin::_time_point_every]

                            # secondary ignore
                            if len(_same_fibers_densities_filtered) < GENERAL_MINIMUM_CORRELATION_TIME_POINTS[EXPERIMENT]:
                                continue

                            _different_correlation = compute_lib.correlation(
                                compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                                compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
                            )

                            _point_distance = compute_lib.distance_from_a_point_to_a_line(
                                _line=[-1, -1, 1, 1],
                                _point=[_same_correlation, _different_correlation]
                            )
                            if _same_correlation > _different_correlation:
                                _y_arrays[_time_point_index].append(_point_distance)
                                _higher_same_counter += 1
                            else:
                                _y_arrays[_time_point_index].append(-_point_distance)

        print('Total points:', len(_y_arrays[_time_point_index]))
        print('Wilcoxon around the zero:')
        print(wilcoxon(_y_arrays[_time_point_index]))
        print('Higher same amount:', _higher_same_counter / len(_y_arrays[_time_point_index]))
        _x_array.append(_time_point_every * TIME_RESOLUTION[EXPERIMENT])

    # plot
    # _colors_array = ['#844b00', '#ea8500', '#edbc80']
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
                'title': 'Distance from y = x',
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
        _filename='plot_' + EXPERIMENT
    )


if __name__ == '__main__':
    main()
