import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths, config
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

TIME_POINT = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STDS = [0, 0.25, 0.5, 0.75]
PAIR_DISTANCES = [5, 7, 9]


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
                'time_points': TIME_POINT
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
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_heterogeneities(_simulations, _stds=STDS)
    _simulations = filtering.by_pair_distances(_simulations, _distances=PAIR_DISTANCES)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations)

    for _std in STDS:
        print('Std:', _std)
        _simulations_std = filtering.by_heterogeneity(_simulations, _std=_std)
        _communicating_correlations = []
        _non_communicating_correlations = []
        _y_arrays = [[] for _i in PAIR_DISTANCES]
        for _pair_distance_index, _pair_distance in enumerate(PAIR_DISTANCES):
            print('Pair distance:', _pair_distance)
            _simulations_std_pair_distance = filtering.by_pair_distance(_simulations_std, _distance=_pair_distance)
            print('Total simulations:', len(_simulations_std_pair_distance))
            _communicating_correlations_array = []
            _non_communicating_correlations_array = []
            for _communicating_index in tqdm(range(len(_simulations_std_pair_distance)), desc='Main loop'):
                _communicating_simulation = _simulations_std_pair_distance[_communicating_index]
                _communicating_left_cell_fiber_densities = _fiber_densities[(_communicating_simulation, 'left_cell')]
                _communicating_right_cell_fiber_densities = _fiber_densities[(_communicating_simulation, 'right_cell')]
                _communicating_correlation = compute_lib.correlation(
                    compute_lib.derivative(_communicating_left_cell_fiber_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_communicating_right_cell_fiber_densities, _n=DERIVATIVE)
                )

                # non-communicating
                _simulations_indices = range(len(_simulations_std_pair_distance))
                for _simulation_index_1 in _simulations_indices:
                    # not the same as the communicating pair
                    if _simulation_index_1 == _communicating_index:
                        continue
                    _simulation_1 = _simulations_std_pair_distance[_simulation_index_1]
                    # after simulation 1 index
                    for _simulation_index_2 in _simulations_indices[_simulation_index_1 + 1:]:
                        # and not the same as the communicating pair
                        if _simulation_index_2 == _communicating_index:
                            continue
                        _simulation_2 = _simulations_std_pair_distance[_simulation_index_2]
                        for _simulation_1_cell_id, _simulation_2_cell_id in product(['left_cell', 'right_cell'],
                                                                                    ['left_cell', 'right_cell']):
                            _simulation_1_fiber_densities = _fiber_densities[(_simulation_1, _simulation_1_cell_id)]
                            _simulation_2_fiber_densities = _fiber_densities[(_simulation_2, _simulation_2_cell_id)]
                            _non_communicating_correlations_array.append(compute_lib.correlation(
                                compute_lib.derivative(_simulation_1_fiber_densities, _n=DERIVATIVE),
                                compute_lib.derivative(_simulation_2_fiber_densities, _n=DERIVATIVE)
                            ))
                            _communicating_correlations_array.append(_communicating_correlation)

            if len(_communicating_correlations_array) > 0:
                print('Total points:', len(_communicating_correlations_array))
                _communicating_minus_non_communicating = \
                    np.array(_communicating_correlations_array) - np.array(_non_communicating_correlations_array)
                print('Wilcoxon of communicating minus non-communicating around the zero:')
                print(wilcoxon(_communicating_minus_non_communicating))
                print('Higher communicating amount:', (_communicating_minus_non_communicating > 0).sum() /
                      len(_communicating_minus_non_communicating))

            _communicating_correlations.append(_communicating_correlations_array)
            _non_communicating_correlations.append(_non_communicating_correlations_array)

        # box plot
        _box_y_arrays = [[] for _i in PAIR_DISTANCES]
        for _communicating_correlations_array, _non_communicating_correlations_array, _pair_distance_index in \
                zip(_communicating_correlations, _non_communicating_correlations, range(len(PAIR_DISTANCES))):
            for _x, _y in zip(_communicating_correlations_array, _non_communicating_correlations_array):
                _point_distance = compute_lib.distance_from_a_point_to_a_line(_line=[-1, -1, 1, 1], _point=[_x, _y])
                if _x > _y:
                    _box_y_arrays[_pair_distance_index].append(_point_distance)
                else:
                    _box_y_arrays[_pair_distance_index].append(-_point_distance)

        _colors_array = config.colors(len(PAIR_DISTANCES))
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_pair_distance,
                    boxpoints=False,
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
                ) for _y, _pair_distance, _color in zip(_box_y_arrays, PAIR_DISTANCES, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Pair distance (cell diameter)',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Communicating minus non-communicating correlation',
                    'zeroline': False,
                    'range': [-1, 1],
                    'tickmode': 'array',
                    'tickvals': [-0.8, -0.4, 0, 0.4, 0.8]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_std_' + str(_std)
        )


if __name__ == '__main__':
    main()
