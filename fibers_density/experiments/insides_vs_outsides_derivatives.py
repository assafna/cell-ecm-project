import os
from itertools import product

import numpy as np
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
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['I', 'I\'', 'I\'\'']


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_real_cells(_experiments)
    if BAND:
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

        for _cell_id, _direction in product(['left_cell', 'right_cell'], ['inside', 'outside']):
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

    _x_arrays = [[] for _i in DERIVATIVES]
    _y_arrays = [[] for _i in DERIVATIVES]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple

        # check if directions exist for both cells
        _missing_information = False
        for _cell_id, _direction in product(['left_cell', 'right_cell'], ['inside', 'outside']):
            if (_experiment, _series_id, _group, _cell_id, _direction) not in _rois_dictionary:
                _missing_information = True
                break
        if _missing_information:
            continue

        _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
        _properties = load.group_properties(_experiment, _series_id, _group)
        _insides_arrays = None
        _outsides_arrays = None
        for _direction in ['inside', 'outside']:
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
                break

            if _direction == 'inside':
                _insides_arrays = [_left_cell_fibers_densities, _right_cell_fibers_densities]
            else:
                _outsides_arrays = [_left_cell_fibers_densities, _right_cell_fibers_densities]

        # enough information
        if _insides_arrays is not None and _outsides_arrays is not None:
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _x_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_insides_arrays[0], _n=_derivative),
                    compute_lib.derivative(_insides_arrays[1], _n=_derivative)
                ))
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_outsides_arrays[0], _n=_derivative),
                    compute_lib.derivative(_outsides_arrays[1], _n=_derivative)
                ))

    print('Total pairs:', len(_x_arrays[0]))
    print('Wilcoxon of insides minus outsides around the zero')
    for _x_array, _y_array, _derivative in zip(_x_arrays, _y_arrays, DERIVATIVES):
        _x_minus_y = np.array(_x_array) - np.array(_y_array)
        print('Derivative:', _derivative, wilcoxon(_x_minus_y))

    # plot
    _colors_array = ['#844b00', '#ea8500', '#edbc80']
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=_y,
                name='Derivative ' + str(_derivative),
                mode='markers',
                marker={
                    'size': 15,
                    'color': _color
                },
                opacity=0.7,
                showlegend=False
            ) for _x, _y, _derivative, _color in zip(_x_arrays, _y_arrays, DERIVATIVES, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Insides correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Outsides correlation',
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
            ],
            'annotations': [
                go.layout.Annotation(
                    x=np.mean(_x),
                    y=np.max(_y) + 0.15,
                    text=_text,
                    showarrow=False,
                    font={
                        'size': 25,
                        'color': 'black'
                    }
                ) for _x, _y, _text in zip(_x_arrays, _y_arrays, DERIVATIVES_TEXT)
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
