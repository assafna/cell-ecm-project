import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

OFFSETS_X = [0, 0.5, 1, 1.5]
OFFSET_Y = 0
CELLS_DISTANCE = 5
SIMULATIONS_TIME_POINTS = 50
DERIVATIVE = 2
STD = 0.5


def compute_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x, _cell_id, _direction in \
                product(OFFSETS_X, ['left_cell', 'right_cell'], ['inside', 'outside']):
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': _offset_x,
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
            _fibers_densities[(_keys['simulation'], _keys['offset_x'], _keys['cell_id'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_fibers_densities(_simulations)

    _x_arrays = [[] for _i in OFFSETS_X]
    _y_arrays = [[] for _i in OFFSETS_X]
    for _offset_x_index, _offset_x in enumerate(OFFSETS_X):
        print('Offset X:', _offset_x)
        for _simulation in tqdm(_simulations, desc='Simulations loop'):
            for _direction in ['inside', 'outside']:
                _left_cell_fibers_densities = _fibers_densities[(_simulation, _offset_x, 'left_cell', _direction)]
                _right_cell_fibers_densities = _fibers_densities[(_simulation, _offset_x, 'right_cell', _direction)]
                _correlation = compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_right_cell_fibers_densities, _n=DERIVATIVE)
                )
                if _direction == 'inside':
                    _x_arrays[_offset_x_index].append(_correlation)
                else:
                    _y_arrays[_offset_x_index].append(_correlation)
        print('Wilcoxon of insides minus outsides around the zero:')
        print(wilcoxon(np.array(_x_arrays[_offset_x_index]) - np.array(_y_arrays[_offset_x_index])))

    # 2d plots
    _colors_array = ['#011f4b', '#00417c', '#2e82bf', '#56caed']
    _legendgroup_array = ['group_1', 'group_1', 'group_2', 'group_2']
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=_y,
                name='Dist. ' + str(_offset_x),
                mode='markers',
                marker={
                    'size': 10,
                    'color': _color
                },
                legendgroup=_legendgroup
            ) for _x, _y, _offset_x, _color, _legendgroup in
            zip(_x_arrays, _y_arrays, OFFSETS_X, _colors_array, _legendgroup_array)
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
            'legend': {
                'x': 0.1,
                'y': 1,
                'xanchor': 'left',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white',
                'orientation': 'h'
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

    for _x, _y, _offset_x, _color in zip(_x_arrays, _y_arrays, OFFSETS_X, _colors_array):
        _fig = go.Figure(
            data=go.Scatter(
                x=_x,
                y=_y,
                name='Dist. ' + str(_offset_x),
                mode='markers',
                marker={
                    'size': 10,
                    'color': _color
                },
                showlegend=True
            ),
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
                'legend': {
                    'xanchor': 'left',
                    'x': 0.1,
                    'yanchor': 'top',
                    'bordercolor': 'black',
                    'borderwidth': 2,
                    'bgcolor': 'white'
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
            _filename='plot_offset_x_' + str(_offset_x)
        )

    # box plot
    _box_y_arrays = [[] for _i in OFFSETS_X]
    for _x_array, _y_array, _offset_x_index in zip(_x_arrays, _y_arrays, range(len(OFFSETS_X))):
        for _x, _y in zip(_x_array, _y_array):
            _point_distance = compute_lib.distance_from_a_point_to_a_line(_line=[-1, -1, 1, 1], _point=[_x, _y])
            if _x > _y:
                _box_y_arrays[_offset_x_index].append(_point_distance)
            else:
                _box_y_arrays[_offset_x_index].append(-_point_distance)

    _fig = go.Figure(
        data=[
            go.Box(
                y=_y,
                name=str(_offset_x),
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
                showlegend=False
            ) for _y, _offset_x, _color in zip(_box_y_arrays, OFFSETS_X, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Distance from cell (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': OFFSETS_X,
                'type': 'category'
            },
            'yaxis': {
                'title': 'Distance from y = x',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-0.2, 0, 0.2, 0.4]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_box'
    )


if __name__ == '__main__':
    main()
