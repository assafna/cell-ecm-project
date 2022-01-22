import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

# by low connectivity
TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
# by low connectivity
STD = {
    False: 0.5,
    True: 0.25
}
PAIR_DISTANCE = [4, 5, 7, 9]


def compute_fiber_densities(_simulations, _low_connectivity):
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
                'time_points': TIME_POINT[_low_connectivity]
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main(_low_connectivity=False):
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT[_low_connectivity])
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD[_low_connectivity])
    _simulations = filtering.by_pair_distances(_simulations, _distances=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations, _low_connectivity)

    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in tqdm(range(len(_simulations)), desc='Main loop'):
        _same_simulation = _simulations[_same_index]
        _same_left_cell_fiber_densities = _fiber_densities[(_same_simulation, 'left_cell')]
        _same_right_cell_fiber_densities = _fiber_densities[(_same_simulation, 'right_cell')]
        _same_correlation = compute_lib.correlation(
            compute_lib.derivative(_same_left_cell_fiber_densities, _n=DERIVATIVE),
            compute_lib.derivative(_same_right_cell_fiber_densities, _n=DERIVATIVE)
        )
        for _different_index in range(len(_simulations)):
            if _same_index != _different_index:
                _different_simulation = _simulations[_different_index]
                for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                 ['left_cell', 'right_cell']):
                    _same_fiber_densities = \
                        _fiber_densities[(_same_simulation, _same_cell_id)]
                    _different_fiber_densities = \
                        _fiber_densities[(_different_simulation, _different_cell_id)]
                    _different_correlations_array.append(compute_lib.correlation(
                        compute_lib.derivative(_same_fiber_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_different_fiber_densities, _n=DERIVATIVE)
                    ))
                    _same_correlations_array.append(_same_correlation)

    print('Total points:', len(_same_correlations_array))
    _same_minus_different = \
        np.array(_same_correlations_array) - np.array(_different_correlations_array)
    print('Wilcoxon of same minus different around the zero:')
    print(wilcoxon(_same_minus_different))
    print('Higher same amount:', (_same_minus_different > 0).sum() /
          len(_same_minus_different))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_same_correlations_array,
            y=_different_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#2e82bf'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Same network correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Different network correlation',
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
