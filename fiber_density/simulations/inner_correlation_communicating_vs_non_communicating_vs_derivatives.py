import os
from itertools import product
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import ranksums
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, config, compute, paths
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCE = 5
TIME_POINTS = 50
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


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
            _fiber_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINTS)
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

    _derivatives_communicating = [[] for _i in DERIVATIVES]
    _derivatives_non_communicating = [[] for _i in DERIVATIVES]
    # communicating loop
    for _simulation in tqdm(_simulations, desc='Communicating pairs loop'):
        _left_cell_fiber_densities = _fiber_densities[(_simulation, 'left_cell')]
        _right_cell_fiber_densities = _fiber_densities[(_simulation, 'right_cell')]

        # derivatives
        for _derivative_index, _derivative in enumerate(DERIVATIVES):
            _correlation = compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities, _n=_derivative),
                compute_lib.derivative(_right_cell_fiber_densities, _n=_derivative)
            )
            _derivatives_communicating[_derivative_index].append(_correlation)

    # non-communicating loop
    _simulations_indices = range(len(_simulations))
    for _simulation_1_index in tqdm(_simulations_indices, desc='Non-communicating pairs loop'):
        _simulation_1 = _simulations[_simulation_1_index]
        for _simulation_2_index in _simulations_indices[_simulation_1_index + 1:]:
            _simulation_2 = _simulations[_simulation_2_index]
            for _simulation_1_cell_id, _simulation_2_cell_id in product(['left_cell', 'right_cell'],
                                                                        ['left_cell', 'right_cell']):
                _simulation_1_fiber_densities = _fiber_densities[(_simulation_1, _simulation_1_cell_id)]
                _simulation_2_fiber_densities = _fiber_densities[(_simulation_2, _simulation_2_cell_id)]

                # derivatives
                for _derivative_index, _derivative in enumerate(DERIVATIVES):
                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_simulation_1_fiber_densities, _n=_derivative),
                        compute_lib.derivative(_simulation_2_fiber_densities, _n=_derivative)
                    )
                    _derivatives_non_communicating[_derivative_index].append(_correlation)

    # rank sums
    print('Wilcoxon rank-sum tests between communicating and non-communicating:')
    for _derivative_index, _derivative in enumerate(DERIVATIVES):
        print('Derivative: ', _derivative, ranksums(
            _derivatives_communicating[_derivative_index], _derivatives_non_communicating[_derivative_index]))

    # plot
    _data = []
    _colors_array = config.colors(2)
    for _communicating, _communicating_text, _derivatives, _color in \
            zip([True, False], ['Communicating', 'Non-communicating'],
                [_derivatives_communicating, _derivatives_non_communicating],
                _colors_array):
        _y = []
        _x = []
        for _derivative_index, _derivative in enumerate(DERIVATIVES_TEXT):
            _y += _derivatives[_derivative_index]
            _x += [_derivative for _i in _derivatives[_derivative_index]]
        _data.append(
            go.Box(
                y=_y,
                x=_x,
                name=_communicating_text,
                boxpoints='all' if _communicating else False,
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
                opacity=0.7
            )
        )

    _fig = go.Figure(
        data=_data,
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
            },
            'boxmode': 'group',
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
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
