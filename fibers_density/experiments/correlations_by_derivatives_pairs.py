import os
from itertools import product

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

EXPERIMENTS = ['SN16']
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
BAND = True
OUT_OF_BOUNDARIES = False
CELLS_DISTANCE_RANGE = [4, 10]
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def main(_directions=None):
    if _directions is None:
        _directions = ['inside', 'outside']

    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_real_cells(_experiments)
    _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
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

        for _cell_id, _direction in product(['left_cell', 'right_cell'], _directions):
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
                'direction': _direction,
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'direction'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    for _direction in _directions:
        _y_arrays = [[] for _i in DERIVATIVES]
        for _tuple in tqdm(_experiments, desc='Experiments loop'):
            _experiment, _series_id, _group = _tuple

            if (_experiment, _series_id, _group, 'left_cell', _direction) not in _rois_dictionary or \
                    (_experiment, _series_id, _group, 'right_cell', _direction) not in _rois_dictionary:
                continue

            _properties = load.group_properties(_experiment, _series_id, _group)

            _left_cell_fibers_densities = \
                _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell', _direction)]
            _right_cell_fibers_densities = \
                _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell', _direction)]

            _left_cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fibers_densities)
            _right_cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fibers_densities)

            if not OUT_OF_BOUNDARIES:
                _left_cell_fibers_densities, _right_cell_fibers_densities = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _left_cell_fibers_densities, _right_cell_fibers_densities
                    )
            else:
                _left_cell_fibers_densities = [_fibers_density[0] for _fibers_density in _left_cell_fibers_densities]
                _right_cell_fibers_densities = [_fibers_density[0] for _fibers_density in _right_cell_fibers_densities]

            # ignore small arrays
            if len(_left_cell_fibers_densities) < MINIMUM_CORRELATION_TIME_POINTS[_experiment] or \
                    len(_right_cell_fibers_densities) < MINIMUM_CORRELATION_TIME_POINTS[_experiment]:
                continue

            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fibers_densities, _n=_derivative),
                    compute_lib.derivative(_right_cell_fibers_densities, _n=_derivative)
                ))

        print('Direction:', _direction)
        print('Total pairs:', len(_y_arrays[0]))
        print('Wilcoxon around the zero')
        for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
            print('Derivative:', _derivative, wilcoxon(_y_array))

        # plot
        _y_title = 'Inner correlation' if _direction == 'inside' else 'Outer correlation'
        _colors_array = ['#844b00', '#ea8500', '#edbc80']
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_derivative,
                    boxpoints='all',
                    jitter=1,
                    pointpos=0,
                    line={
                        'width': 1
                    },
                    fillcolor='white',
                    marker={
                        'size': 10,
                        'color': _color
                    },
                    opacity=0.7,
                    showlegend=False
                ) for _y, _derivative, _color in zip(_y_arrays, DERIVATIVES_TEXT, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Density by derivatives',
                    'zeroline': False
                },
                'yaxis': {
                    'title': _y_title,
                    'range': [-1, 1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_direction_' + _direction
        )


if __name__ == '__main__':
    main()
