import os
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCE = 7
TIME_POINTS = 50
TIME_POINTS_STEP = 3

DERIVATIVES = [0, 1, 2]
DERIVATIVES_COLORS = {
    0: '#011f4b',
    1: '#005b96',
    2: '#74c2e8'
}
DERIVATIVES_Y_TICKVALS = {
    0: [0, 2, 4],
    1: [0, 0.1, 0.2],
    2: [-0.04, 0, 0.04]
}


def compute_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
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


def main():
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

    _y_arrays = [[] for _i in DERIVATIVES]
    for _simulation in tqdm(_simulations, desc='Simulations loop'):
        _normalization = load.normalization(_simulation)
        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fiber_densities = _fiber_densities[(_simulation, _cell_id)]
            _cell_fiber_densities_normalized = compute_lib.z_score_array(
                _array=_cell_fiber_densities,
                _average=_normalization['average'],
                _std=_normalization['std']
            )
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(
                    compute_lib.derivative(_cell_fiber_densities_normalized, _n=_derivative)
                )

    print('Total cells:', len(_y_arrays[0]))

    # plot
    for _derivative, _y_array in zip(DERIVATIVES, _y_arrays):
        _fig = go.Figure(
            data=go.Scatter(
                x=list(range(TIME_POINTS))[::TIME_POINTS_STEP],
                y=np.mean(_y_array, axis=0)[::TIME_POINTS_STEP],
                name='Fiber density (z-score)',
                error_y={
                    'type': 'data',
                    'array': np.std(_y_array, axis=0)[::TIME_POINTS_STEP],
                    'thickness': 1,
                    'color': DERIVATIVES_COLORS[_derivative]
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': DERIVATIVES_COLORS[_derivative]
                },
                opacity=0.7
            ),
            layout={
                'xaxis': {
                    'title': 'Cell contraction (%)',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Fiber density (z-score)',
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': DERIVATIVES_Y_TICKVALS[_derivative]
                },
                'shapes': [
                    {
                        'type': 'line',
                        'x0': -2,
                        'y0': (np.mean(_y_array, axis=0)[0] - np.std(_y_array, axis=0)[0]) * 1.5,
                        'x1': 53,
                        'y1': (np.mean(_y_array, axis=0)[0] - np.std(_y_array, axis=0)[0]) * 1.5,
                        'line': {
                            'color': 'black',
                            'width': 2
                        }
                    },
                    {
                        'type': 'line',
                        'x0': -2,
                        'y0': (np.mean(_y_array, axis=0)[0] - np.std(_y_array, axis=0)[0]) * 1.5,
                        'x1': -2,
                        'y1': (np.mean(_y_array, axis=0)[-1] + np.std(_y_array, axis=0)[-1]),
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
            _filename='plot_derivative_' + str(_derivative)
        )


if __name__ == '__main__':
    main()
