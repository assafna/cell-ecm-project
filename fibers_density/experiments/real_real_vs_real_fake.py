import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}


def main(_offset_y=0.5, _high_time_resolution=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_band(_experiments)
    _tuples_matched = organize.by_matched_real_and_fake(_experiments)
    print('Total matched pairs:', len(_tuples_matched))

    _arguments = []
    for _matched_tuple in _tuples_matched:
        for _tuple in _matched_tuple:
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
                    'offset_y': _offset_y,
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

    _real_real_correlations_array = []
    _real_fake_correlations_array = []
    _valid_real_real_tuples = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]
        _tuples_matched = organize.by_matched_real_and_fake(_experiment_tuples)
        print('Matched pairs:', len(_tuples_matched))

        for _tuple_matched in tqdm(_tuples_matched, desc='Main loop'):
            _real_real_tuple, _real_fake_tuple = _tuple_matched

            _real_real_experiment, _real_real_series, _real_real_group = _real_real_tuple
            _, _, _real_fake_group = _real_fake_tuple

            _real_left_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_real_real_experiment, _real_real_series, _real_real_group, 'left_cell')
                ]
            _real_right_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_real_real_experiment, _real_real_series, _real_real_group, 'right_cell')
                ]
            _fake_left_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_real_real_experiment, _real_real_series, _real_fake_group, 'left_cell')
                ]
            _fake_right_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_real_real_experiment, _real_real_series, _real_fake_group, 'right_cell')
                ]

            _properties = \
                load.group_properties(_real_real_experiment, _real_real_series, _real_real_group)

            _real_left_cell_fibers_densities = compute.remove_blacklist(
                _real_real_experiment,
                _real_real_series,
                _properties['cells_ids']['left_cell'],
                _real_left_cell_fibers_densities
            )
            _real_right_cell_fibers_densities = compute.remove_blacklist(
                _real_real_experiment,
                _real_real_series,
                _properties['cells_ids']['right_cell'],
                _real_right_cell_fibers_densities
            )
            _fake_left_cell_fibers_densities = compute.remove_blacklist(
                _real_real_experiment,
                _real_real_series,
                _properties['cells_ids']['left_cell'],
                _fake_left_cell_fibers_densities
            )
            _fake_right_cell_fibers_densities = compute.remove_blacklist(
                _real_real_experiment,
                _real_real_series,
                _properties['cells_ids']['right_cell'],
                _fake_right_cell_fibers_densities
            )

            _real_left_cell_fibers_densities_filtered, _real_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _real_left_cell_fibers_densities, _real_right_cell_fibers_densities
                )

            # ignore small arrays
            if len(_real_left_cell_fibers_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_real_real_experiment]:
                continue

            _real_real_correlation = compute_lib.correlation(
                compute_lib.derivative(_real_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_real_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )

            for _real_cell_fiber_densities, _fake_cell_fiber_densities in \
                    zip([_real_left_cell_fibers_densities, _real_right_cell_fibers_densities],
                        [_fake_left_cell_fibers_densities, _fake_right_cell_fibers_densities]):
                _real_cell_fiber_densities_filtered, _fake_cell_fiber_densities_filtered = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _real_cell_fiber_densities, _fake_cell_fiber_densities
                    )

                # ignore small arrays
                if len(_real_cell_fiber_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_real_real_experiment]:
                    continue

                _real_fake_correlation = compute_lib.correlation(
                    compute_lib.derivative(_real_cell_fiber_densities_filtered, _n=DERIVATIVE),
                    compute_lib.derivative(_fake_cell_fiber_densities_filtered, _n=DERIVATIVE)
                )

                _real_real_correlations_array.append(_real_real_correlation)
                _real_fake_correlations_array.append(_real_fake_correlation)

                if _real_real_tuple not in _valid_real_real_tuples:
                    _valid_real_real_tuples.append(_real_real_tuple)

    print('Total real-real pairs:', len(_valid_real_real_tuples))
    _real_minus_fake = np.array(_real_real_correlations_array) - np.array(_real_fake_correlations_array)
    print('Wilcoxon of real-real minus real-fake around the zero:')
    print(wilcoxon(_real_minus_fake))
    print('Higher real-real amount:', (_real_minus_fake > 0).sum() / len(_real_minus_fake))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_real_real_correlations_array,
            y=_real_fake_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Real-real correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Real-fake correlation',
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
