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

OFFSET_X = 0
OFFSET_Y = 0
CELLS_DISTANCES = [4, 5, 7, 9]
SIMULATIONS_TIME_POINTS = 50
DERIVATIVE = 2
STD = 0.5


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
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distances(_simulations, _distances=CELLS_DISTANCES)
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_fibers_densities(_simulations)

    _x_arrays = [[] for _i in CELLS_DISTANCES]
    _y_arrays = [[] for _i in CELLS_DISTANCES]
    for _distance_index, _distance in enumerate(CELLS_DISTANCES):
        _distance_simulations = filtering.by_distance(_simulations, _distance=_distance)
        print('Distance:', _distance, 'Total simulations:', len(_distance_simulations))
        for _simulation in tqdm(_distance_simulations, desc='Simulations loop'):
            for _direction in ['inside', 'outside']:
                _left_cell_fibers_densities = _fibers_densities[(_simulation, 'left_cell', _direction)]
                _right_cell_fibers_densities = _fibers_densities[(_simulation, 'right_cell', _direction)]
                _correlation = compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_right_cell_fibers_densities, _n=DERIVATIVE)
                )
                if _direction == 'inside':
                    _x_arrays[_distance_index].append(_correlation)
                else:
                    _y_arrays[_distance_index].append(_correlation)
        print('Wilcoxon of insides minus outsides around the zero:')
        print(wilcoxon(np.array(_x_arrays[_distance_index]) - np.array(_y_arrays[_distance_index])))

    # 2d plots
    _colors_array = ['#011f4b', '#00417c', '#2e82bf', '#56caed']
    _legendgroup_array = ['group_1', 'group_1', 'group_2', 'group_2']
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=_y,
                name='Dist. ' + str(_distance),
                mode='markers',
                marker={
                    'size': 10,
                    'color': _color
                },
                legendgroup=_legendgroup
            ) for _x, _y, _distance, _color, _legendgroup in
            zip(_x_arrays, _y_arrays, CELLS_DISTANCES, _colors_array, _legendgroup_array)
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

    for _x, _y, _distance, _color in zip(_x_arrays, _y_arrays, CELLS_DISTANCES, _colors_array):
        _fig = go.Figure(
            data=go.Scatter(
                x=_x,
                y=_y,
                name='Dist. ' + str(_distance),
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
            _filename='plot_distance_' + str(_distance)
        )

    # box plot
    _box_y_arrays = [[] for _i in CELLS_DISTANCES]
    for _x_array, _y_array, _distance_index in zip(_x_arrays, _y_arrays, range(len(CELLS_DISTANCES))):
        for _x, _y in zip(_x_array, _y_array):
            _point_distance = compute_lib.distance_from_a_point_to_a_line(_line=[-1, -1, 1, 1], _point=[_x, _y])
            if _x > _y:
                _box_y_arrays[_distance_index].append(_point_distance)
            else:
                _box_y_arrays[_distance_index].append(-_point_distance)

    _fig = go.Figure(
        data=[
            go.Box(
                y=_y,
                name=str(_distance),
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
            ) for _y, _distance, _color in zip(_box_y_arrays, CELLS_DISTANCES, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Cells distance (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': CELLS_DISTANCES,
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
