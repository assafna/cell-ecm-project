import os
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, paths
from libs.simulations import config
from libs.simulations import filtering
from libs.simulations import load
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCE = 7
TIME_POINTS = 50
TIME_STEP = 3


def compute_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def compute_data():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_pair_distance(_simulations, _distance=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations)

    _simulations_fiber_densities = [[] for _i in range(TIME_POINTS)]
    for _simulation in tqdm(_simulations, desc='Main loop'):
        _normalization = load.normalization(_simulation)

        for _time_point in range(TIME_POINTS):
            for _cell_id in ['left_cell', 'right_cell']:
                _fiber_density = _fiber_densities[(_simulation, _cell_id)][_time_point]

                _normalized_fiber_density = compute_lib.z_score(
                    _fiber_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _simulations_fiber_densities[_time_point].append(_normalized_fiber_density)

    print('Total pairs:', len(_simulations_fiber_densities[0]))

    return _simulations_fiber_densities


def main():
    _fiber_densities = compute_data()

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=list(range(TIME_POINTS))[::TIME_STEP],
                y=[np.mean(_array) for _array in _fiber_densities][::TIME_STEP],
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _fiber_densities][::TIME_STEP],
                    'thickness': 1,
                    'color': '#005b96'
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#005b96'
                },
                opacity=0.7,
                showlegend=False
            )
        ],
        layout={
            'xaxis': {
                'title': 'Cell contraction (%)',
                'titlefont': {
                    'color': '#005b96'
                },
                'tickfont': {
                    'color': '#005b96'
                },
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 25, 50]
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-1.7, 6],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 2, 4, 6]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -2,
                    'y0': -1.5,
                    'x1': 53,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -2,
                    'y0': -1.5,
                    'x1': -2,
                    'y1': 6,
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
