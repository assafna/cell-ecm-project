import os

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

EXPERIMENTS = ['SN16']
EXPERIMENTS_BLEB = ['SN20_Bleb_fromStart']
TIME_POINT = 18
CELLS_DISTANCE_RANGE = [4, 10]
OFFSET_X_STEP = 0.2
OFFSET_Z = 0
OFFSET_Y = 0
OUT_OF_BOUNDARIES = False


def compute_fibers(_experiments):
    _experiments = load.experiments_groups_as_tuples(_experiments)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_real_cells(_experiments)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    print('Total experiments:', len(_experiments))

    _max_offsets_x = []
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _cells_distance = compute.cells_distance_in_cell_size_time_point(
            _experiment, _series_id, _group, TIME_POINT - 1
        )
        _offsets_x = \
            np.arange(start=0, stop=_cells_distance / 2 - 0.5 - ROI_LENGTH, step=OFFSET_X_STEP)
        if len(_offsets_x) > len(_max_offsets_x):
            _max_offsets_x = _offsets_x
        for _offset_x in _offsets_x:
            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': ROI_LENGTH,
                    'length_y': ROI_HEIGHT,
                    'length_z': ROI_WIDTH,
                    'offset_x': _offset_x,
                    'offset_y': OFFSET_Y,
                    'offset_z': OFFSET_Z,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': TIME_POINT - 1
                })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = [[] for _i in range(len(_max_offsets_x))]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple
        for _offset_x_index, _offset_x in enumerate(_max_offsets_x):
            for _cell_id in ['left_cell', 'right_cell']:
                if (_experiment, _series_id, _group, _offset_x, _cell_id) in _rois_dictionary:
                    _normalization = \
                        load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
                    _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _offset_x, _cell_id)][0]
                    _fibers_density = _fibers_densities[_roi_tuple]

                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    if not np.isnan(_normalized_fibers_density):
                        _experiments_fibers_densities[_offset_x_index].append(_normalized_fibers_density)

    return _experiments_fibers_densities, _max_offsets_x


def main():
    print('Regular experiments')
    _regular_experiments, _regular_offsets_x = compute_fibers(EXPERIMENTS)

    print('Bleb experiments')
    _bleb_experiments, _bleb_offsets_x = compute_fibers(EXPERIMENTS_BLEB)

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_regular_offsets_x,
                y=[np.mean(_array) for _array in _regular_experiments],
                name='No bleb',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _regular_experiments],
                    'thickness': 1,
                    'color': '#005b96'
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#005b96'
                },
                opacity=0.7
            ),
            go.Scatter(
                x=_bleb_offsets_x,
                y=[np.mean(_array) for _array in _bleb_experiments],
                name='Bleb',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _bleb_experiments],
                    'thickness': 1,
                    'color': '#ea8500'
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#ea8500'
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-1.7, 13],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': 3.4,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': -0.2,
                    'y1': 13,
                    'line': {
                        'color': 'black',
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
