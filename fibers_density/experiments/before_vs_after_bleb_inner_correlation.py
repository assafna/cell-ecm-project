import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
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

    _n_pairs = 0
    _before_correlations = []
    _after_correlations = []
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple

        _left_cell_fibers_densities = \
            _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fibers_densities = \
            _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fibers_densities)
        _right_cell_fibers_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fibers_densities)

        _before_left_cell_fiber_densities = \
            _left_cell_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]
        _before_right_cell_fiber_densities = \
            _right_cell_fibers_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_POINT]

        _after_left_cell_fiber_densities = \
            _left_cell_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]
        _after_right_cell_fiber_densities = \
            _right_cell_fibers_densities[AFTER_BLEB_INJECTION_FIRST_TIME_POINT:]

        _before_left_cell_fibers_densities_filtered, _before_right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _before_left_cell_fiber_densities, _before_right_cell_fiber_densities)

        _after_left_cell_fibers_densities_filtered, _after_right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _after_left_cell_fiber_densities, _after_right_cell_fiber_densities)

        # ignore small arrays
        if len(_before_left_cell_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION or \
                len(_after_left_cell_fibers_densities_filtered) < MINIMUM_TIME_POINTS_CORRELATION:
            continue

        _n_pairs += 1

        _before_correlations.append(compute_lib.correlation(
            compute_lib.derivative(_before_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_before_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
        ))
        _after_correlations.append(compute_lib.correlation(
            compute_lib.derivative(_after_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_after_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
        ))

    print('Total pairs:', _n_pairs)
    _before_minus_after = np.array(_before_correlations) - np.array(_after_correlations)
    print('Wilcoxon of before minus after around the zero:')
    print(wilcoxon(_before_minus_after))
    print('Higher before amount:', (_before_minus_after > 0).sum() / len(_before_minus_after))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_before_correlations,
            y=_after_correlations,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Correlation before bleb',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Correlation after bleb',
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
