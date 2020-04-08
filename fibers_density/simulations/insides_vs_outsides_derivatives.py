import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
CELLS_DISTANCE = 7
SIMULATIONS_TIME_POINTS = 50
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['I', 'I\'', 'I\'\'']


def compute_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id, _direction in product(['left_cell', 'right_cell'], ['inside', 'outside']):
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': _direction,
                'time_points': SIMULATIONS_TIME_POINTS
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[(_keys['simulation'], _keys['cell_id'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_fibers_densities(_simulations)

    _x_arrays = [[] for _i in DERIVATIVES]
    _y_arrays = [[] for _i in DERIVATIVES]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _normalization = load.normalization(_simulation)
        _normalization = [_normalization['average'], _normalization['std']]
        for _direction in ['inside', 'outside']:
            _left_cell_fibers_densities = _fibers_densities[(_simulation, 'left_cell', _direction)]
            _right_cell_fibers_densities = _fibers_densities[(_simulation, 'right_cell', _direction)]
            _left_cell_fibers_densities_normalized = \
                compute_lib.z_score_fibers_densities_array(_left_cell_fibers_densities, _normalization)
            _right_cell_fibers_densities_normalized = \
                compute_lib.z_score_fibers_densities_array(_right_cell_fibers_densities, _normalization)
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _correlation = compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fibers_densities_normalized, _n=_derivative),
                    compute_lib.derivative(_right_cell_fibers_densities_normalized, _n=_derivative)
                )
                if _direction == 'inside':
                    _x_arrays[_derivative_index].append(_correlation)
                else:
                    _y_arrays[_derivative_index].append(_correlation)

    # plot
    _colors_array = ['#011f4b', '#005b96', '#74c2e8']
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
