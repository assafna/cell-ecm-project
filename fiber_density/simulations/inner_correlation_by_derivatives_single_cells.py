import os

import numpy as np
import plotly.graph_objs as go
from multiprocess.pool import Pool
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, config, compute, paths
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
TIME_POINTS = 50
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_simulations_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _direction in ['left', 'right', 'up', 'down']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                'offset_y': OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_simulations_fiber_densities(_simulations)

    _y_arrays = [[] for _i in DERIVATIVES]
    for _index_1 in tqdm(range(len(_simulations)), desc='Simulations loop'):
        _simulation_1 = _simulations[_index_1]
        _cell_1_fiber_densities = \
            [_fiber_densities[(_simulation_1, _direction)] for _direction in ['left', 'right', 'up', 'down']]
        _cell_1_fiber_densities = np.mean(_cell_1_fiber_densities, axis=0)
        for _index_2 in range(_index_1 + 1, len(_simulations)):
            _simulation_2 = _simulations[_index_2]
            _cell_2_fiber_densities = \
                [_fiber_densities[(_simulation_2, _direction)] for _direction in ['left', 'right', 'up', 'down']]
            _cell_2_fiber_densities = np.mean(_cell_2_fiber_densities, axis=0)
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_cell_1_fiber_densities, _n=_derivative),
                    compute_lib.derivative(_cell_2_fiber_densities, _n=_derivative)
                ))

    print('Total points:', len(_y_arrays[0]))
    print('Wilcoxon around the zero')
    for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
        print('Derivative:', _derivative, wilcoxon(_y_array))

    # plot
    _colors_array = config.colors(3)
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
                'title': 'Fiber density derivative',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Correlation',
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
