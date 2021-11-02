import os
from itertools import product
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import ranksums
from tqdm import tqdm

from libs import config_lib, compute_lib
from libs.simulations import config, compute, load, filtering, paths
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCES = [4, 5, 7, 9]
TIME_POINTS = 50
DERIVATIVE = 2
STD = 0.5


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
    with Pool(config_lib.CPUS_TO_USE) as _p:
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
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_pair_distances(_simulations, _distances=PAIR_DISTANCES)
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations)

    _pair_distances_communicating = [[] for _i in PAIR_DISTANCES]
    _pair_distances_non_communicating = [[] for _i in PAIR_DISTANCES]
    # pair distances loop
    for _pair_distance_index, _pair_distance in enumerate(PAIR_DISTANCES):
        print('Pair distance:', _pair_distance)
        _pair_distance_simulations = filtering.by_pair_distance(_simulations, _distance=_pair_distance)

        # communicating loop
        for _simulation in tqdm(_pair_distance_simulations, desc='Communicating loop'):
            _left_cell_fiber_densities = _fiber_densities[(_simulation, 'left_cell')]
            _right_cell_fiber_densities = _fiber_densities[(_simulation, 'right_cell')]
            _correlation = compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fiber_densities, _n=DERIVATIVE)
            )
            _pair_distances_communicating[_pair_distance_index].append(_correlation)

        # non-communicating loop
        _simulations_indices = range(len(_pair_distance_simulations))
        for _simulation_1_index in tqdm(_simulations_indices, desc='Non-communicating pairs loop'):
            _simulation_1 = _pair_distance_simulations[_simulation_1_index]
            for _simulation_2_index in _simulations_indices[_simulation_1_index + 1:]:
                _simulation_2 = _pair_distance_simulations[_simulation_2_index]
                for _simulation_1_cell_id, _simulation_2_cell_id in product(['left_cell', 'right_cell'],
                                                                            ['left_cell', 'right_cell']):
                    _simulation_1_fiber_densities = _fiber_densities[(_simulation_1, _simulation_1_cell_id)]
                    _simulation_2_fiber_densities = _fiber_densities[(_simulation_2, _simulation_2_cell_id)]
                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_simulation_1_fiber_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_simulation_2_fiber_densities, _n=DERIVATIVE)
                    )
                    _pair_distances_non_communicating[_pair_distance_index].append(_correlation)

        # rank sums
        print('Wilcoxon rank-sum tests between communicating and non-communicating:',
              ranksums(_pair_distances_communicating[_pair_distance_index],
                       _pair_distances_non_communicating[_pair_distance_index]))

    # plot
    _data = []
    _colors_array = config.colors(2)
    for _communicating, _communicating_text, _pair_distances, _color in \
            zip([True, False], ['Communicating', 'Non-communicating'],
                [_pair_distances_communicating, _pair_distances_non_communicating],
                _colors_array):
        _y = []
        _x = []
        for _pair_distance_index, _pair_distance in enumerate(PAIR_DISTANCES):
            _y += _pair_distances[_pair_distance_index]
            _x += [_pair_distance for _i in _pair_distances[_pair_distance_index]]
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
                'title': 'Pair distance (cell diameter)',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': PAIR_DISTANCES,
                'type': 'category'
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
