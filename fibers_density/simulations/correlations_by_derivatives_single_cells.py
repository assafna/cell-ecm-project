import os

import numpy as np
import plotly.graph_objs as go
from multiprocess.pool import Pool
from scipy.stats import wilcoxon
from statsmodels.tsa.stattools import kpss, adfuller
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


def compute_simulations_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _direction in ['left', 'right', 'up', 'down']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': config.ROI_HEIGHT
                if _direction in ['up', 'down'] else config.ROI_WIDTH,
                'length_y': config.ROI_WIDTH
                if _direction in ['up', 'down'] else config.ROI_HEIGHT,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                'offset_y': OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_POINTS
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[(_keys['simulation'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    print('Total simulations:', len(_simulations))

    _fibers_densities = compute_simulations_fibers_densities(_simulations)

    _y_arrays = [[] for _i in DERIVATIVES]
    for _simulation_1_index in tqdm(range(len(_simulations)), desc='Simulations loop'):
        _simulation_1 = _simulations[_simulation_1_index]
        _cell_1_fibers_densities = \
            [_fibers_densities[(_simulation_1, _direction)] for _direction in ['left', 'right', 'up', 'down']]
        _cell_1_fibers_densities = np.mean(_cell_1_fibers_densities, axis=0)
        for _simulation_2_index in range(_simulation_1_index + 1, len(_simulations)):
            _simulation_2 = _simulations[_simulation_2_index]
            _cell_2_fibers_densities = \
                [_fibers_densities[(_simulation_2, _direction)] for _direction in ['left', 'right', 'up', 'down']]
            _cell_2_fibers_densities = np.mean(_cell_2_fibers_densities, axis=0)
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_cell_1_fibers_densities, _n=_derivative),
                    compute_lib.derivative(_cell_2_fibers_densities, _n=_derivative)
                ))

    print('Total points:', len(_y_arrays[0]))
    print('Wilcoxon around the zero')
    for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
        print('Derivative:', _derivative, wilcoxon(_y_array))

    # plot
    _colors_array = ['#011f4b', '#005b96', '#74c2e8']
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
