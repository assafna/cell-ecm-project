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
PAIR_DISTANCE = [5, 7, 9]


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

    _communicating_correlations_array = []
    _non_communicating_correlations_array = []
    for _communicating_index in tqdm(range(len(_simulations)), desc='Main loop'):
        _communicating_simulation = _simulations[_communicating_index]
        _communicating_left_cell_fiber_densities = _fiber_densities[(_communicating_simulation, 'left_cell')]
        _communicating_right_cell_fiber_densities = _fiber_densities[(_communicating_simulation, 'right_cell')]
        _communicating_correlation = compute_lib.correlation(
            compute_lib.derivative(_communicating_left_cell_fiber_densities, _n=DERIVATIVE),
            compute_lib.derivative(_communicating_right_cell_fiber_densities, _n=DERIVATIVE)
        )

        # non-communicating
        _simulations_indices = range(len(_simulations))
        for _simulation_index_1 in _simulations_indices:
            # not the same as the communicating pair
            if _simulation_index_1 == _communicating_index:
                continue
            _simulation_1 = _simulations[_simulation_index_1]
            # after simulation 1 index
            for _simulation_index_2 in _simulations_indices[_simulation_index_1 + 1:]:
                # and not the same as the communicating pair
                if _simulation_index_2 == _communicating_index:
                    continue
                _simulation_2 = _simulations[_simulation_index_2]
                for _simulation_1_cell_id, _simulation_2_cell_id in product(['left_cell', 'right_cell'],
                                                                            ['left_cell', 'right_cell']):
                    _simulation_1_fiber_densities = _fiber_densities[(_simulation_1, _simulation_1_cell_id)]
                    _simulation_2_fiber_densities = _fiber_densities[(_simulation_2, _simulation_2_cell_id)]
                    _non_communicating_correlations_array.append(compute_lib.correlation(
                        compute_lib.derivative(_simulation_1_fiber_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_simulation_2_fiber_densities, _n=DERIVATIVE)
                    ))
                    _communicating_correlations_array.append(_communicating_correlation)

    print('Total points:', len(_communicating_correlations_array))
    _communicating_minus_non_communicating = \
        np.array(_communicating_correlations_array) - np.array(_non_communicating_correlations_array)
    print('Wilcoxon of communicating minus non-communicating around the zero:')
    print(wilcoxon(_communicating_minus_non_communicating))
    print('Higher communicating amount:', (_communicating_minus_non_communicating > 0).sum() /
          len(_communicating_minus_non_communicating))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_communicating_correlations_array,
            y=_non_communicating_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#2e82bf'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Communicating correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Non-communicating correlation',
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
