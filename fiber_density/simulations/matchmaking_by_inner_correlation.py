import os
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

TIME_POINT = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
PAIR_DISTANCE = 4
STD = 0.5
MAX_RANK = 7


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
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_pair_distance(_simulations, _distance=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations)

    _n = 0
    _cells_potential_matches = []
    _cells_ranks = []
    for _simulation_1 in tqdm(_simulations, desc='Main loop'):
        for _cell_1_id, _cell_1_correct_match_cell_id in zip(['left_cell', 'right_cell'], ['right_cell', 'left_cell']):
            _cell_1 = (_simulation_1, _cell_1_id)
            _cell_1_correct_match = (_simulation_1, _cell_1_correct_match_cell_id)
            _cell_1_fiber_densities = _fiber_densities[_cell_1]

            _cell_1_correlations = []
            _cell_1_correct_match_correlation = None
            for _simulation_2 in _simulations:
                for _cell_2_id in ['left_cell', 'right_cell']:
                    _cell_2 = (_simulation_2, _cell_2_id)

                    # same cell
                    if _cell_1 == _cell_2:
                        continue

                    _cell_2_fiber_densities = _fiber_densities[_cell_2]

                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_cell_1_fiber_densities, _n=DERIVATIVE),
                        compute_lib.derivative(_cell_2_fiber_densities, _n=DERIVATIVE)
                    )

                    _cell_1_correlations.append(_correlation)

                    # correct match
                    if _cell_2 == _cell_1_correct_match:
                        _cell_1_correct_match_correlation = _correlation

            # correct match does not exist
            if _cell_1_correct_match_correlation is None:
                continue

            # check matchmaking
            _cell_1_total_potential_matches = len(_cell_1_correlations)
            if _cell_1_total_potential_matches > 1:
                _cell_1_correct_match_rank = 1
                for _potential_match_correlation in sorted(_cell_1_correlations, reverse=True):
                    if _cell_1_correct_match_correlation == _potential_match_correlation:
                        break
                    _cell_1_correct_match_rank += 1

                _n += 1
                _cells_potential_matches.append(_cell_1_total_potential_matches)
                _cells_ranks.append(_cell_1_correct_match_rank)

    # results
    _mean_cells_potential_matches = float(np.mean(_cells_potential_matches))
    _mean_correct_match_probability = 1 / _mean_cells_potential_matches
    _first_place_correct_matches = sum([1 for _rank in _cells_ranks if _rank == 1])
    _first_place_fraction = _first_place_correct_matches / _n

    print('Matchmaking results:')
    print('Total cells:', _n)
    print('Average potential matches per cell:', round(_mean_cells_potential_matches, 2))
    print('Average correct match probability:', round(_mean_correct_match_probability, 2))
    print('Fraction of first place correct matches:', round(_first_place_fraction, 2))

    # plot
    _x = list(range(MAX_RANK))
    _x_text = [str(_rank + 1) for _rank in _x[:-1]] + [str(MAX_RANK) + '+']
    _ranks_sums = [0 for _rank in _x]
    for _rank in _cells_ranks:
        if _rank < MAX_RANK:
            _ranks_sums[_rank - 1] += 1
        else:
            _ranks_sums[-1] += 1
    _y = np.array(_ranks_sums) / _n

    _fig = go.Figure(
        data=go.Bar(
            x=_x,
            y=_y,
            marker={
                'color': '#2e82bf'
            }
        ),
        layout={
            'xaxis': {
                'title': 'Correct match correlation rank',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': _x,
                'ticktext': _x_text
            },
            'yaxis': {
                'title': 'Fraction',
                'range': [0, 1.1],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )

    # correct match probability plot
    _y = [_mean_correct_match_probability] * (MAX_RANK - 1) + [_mean_correct_match_probability * MAX_RANK]
    _fig = go.Figure(
        data=go.Bar(
            x=_x,
            y=_y,
            marker={
                'color': '#2e82bf'
            }
        ),
        layout={
            'xaxis': {
                'title': 'Correct match correlation rank',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': _x,
                'ticktext': _x_text
            },
            'yaxis': {
                'title': 'Fraction',
                'range': [0, 1.1],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_correct_match_prob'
    )


if __name__ == '__main__':
    main()
