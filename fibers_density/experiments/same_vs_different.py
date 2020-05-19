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
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
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
    'SN44': 50,
    'SN45': 50
}


def compute_fibers_densities(_band=True, _high_time_resolution=False, _cells_distance_range=None):
    if _cells_distance_range is None:
        _cells_distance_range = CELLS_DISTANCE_RANGE

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, _cells_distance_range)
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

    _tuples_by_experiment = organize.by_experiment(_experiments)

    _same_correlations_array = []
    _different_correlations_array = []
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

            _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fibers_densities, _same_right_cell_fibers_densities
                )

            # ignore small arrays
            if len(_same_left_cell_fibers_densities_filtered) < \
                    MINIMUM_CORRELATION_TIME_POINTS[_same_experiment]:
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
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

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
                        )

                        _same_correlations_array.append(_same_correlation)
                        _different_correlations_array.append(_different_correlation)

    print('Total points:', len(_same_correlations_array))
    _same_minus_different = \
        np.array(_same_correlations_array) - np.array(_different_correlations_array)
    print('Wilcoxon of same minus different around the zero:')
    print(wilcoxon(_same_minus_different))
    print('Higher same amount:', (_same_minus_different > 0).sum() /
          len(_same_minus_different))

    return _same_correlations_array, _different_correlations_array


def main(_band=True, _high_time_resolution=False, _cells_distance_range=None):
    _same_correlations_array, _different_correlations_array =\
        compute_fibers_densities(_band, _high_time_resolution, _cells_distance_range)

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_same_correlations_array,
            y=_different_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Same network correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Different network correlation',
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
        _filename='plot_band_' + str(_band) + '_high_time_res_' + str(_high_time_resolution) + '_range_' +
                  '_'.join([str(_distance) for _distance in _cells_distance_range])
    )


if __name__ == '__main__':
    main()
