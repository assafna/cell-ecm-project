import os

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
PAIR_DISTANCE_RANGE = [4, 10]
STATIC = False
MINIMUM_CORRELATION_TIME_FRAMES = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}
MAX_RANK = 7


def main(_real_cells=True, _offset_y=0.5, _high_time_resolution=False, _band=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_real_pairs(_experiments, _real_pairs=_real_cells)
    _experiments = filtering.by_fake_static_pairs(_experiments, _fake_static_pairs=STATIC)
    _experiments = filtering.by_band(_experiments, _band=_band)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_frame = len(_properties['time_points'])
        for _time_frame in range(len(_properties['time_points'])):
            _pair_distance = \
                compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame)
            if _pair_distance - 1 - OFFSET_X * 2 < QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER * 2:
                _latest_time_frame = _time_frame - 1
                break

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _tuples_by_experiment = organize.by_experiment(_experiments)

    _n = 0
    _cells_potential_matches = []
    _cells_ranks = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _tuple_1 in tqdm(_experiment_tuples, desc='Main loop'):
            _experiment_1, _series_id_1, _group_1 = _tuple_1
            _experiment_1_properties = load.group_properties(_experiment_1, _series_id_1, _group_1)

            for _cell_1_id, _cell_1_correct_match_cell_id in \
                    zip(['left_cell', 'right_cell'], ['right_cell', 'left_cell']):
                _cell_1 = (_experiment_1, _series_id_1, _group_1, _cell_1_id)
                _cell_1_correct_match = (_experiment_1, _series_id_1, _group_1, _cell_1_correct_match_cell_id)
                _cell_1_fiber_densities = \
                    _experiments_fiber_densities[(_experiment_1, _series_id_1, _group_1, _cell_1_id)]
                _cell_1_fiber_densities = compute.remove_blacklist(
                    _experiment_1,
                    _series_id_1,
                    _experiment_1_properties['cells_ids'][_cell_1_id],
                    _cell_1_fiber_densities
                )

                _cell_1_correlations = []
                _cell_1_correct_match_correlation = None
                for _tuple_2 in _experiment_tuples:
                    _experiment_2, _series_id_2, _group_2 = _tuple_2
                    _experiment_2_properties = load.group_properties(_experiment_2, _series_id_2, _group_2)

                    for _cell_2_id in ['left_cell', 'right_cell']:
                        _cell_2 = (_experiment_2, _series_id_2, _group_2, _cell_2_id)

                        # same cell
                        if _cell_1 == _cell_2:
                            continue

                        _cell_2_fiber_densities = \
                            _experiments_fiber_densities[(_experiment_2, _series_id_2, _group_2, _cell_2_id)]
                        _cell_2_fiber_densities = compute.remove_blacklist(
                            _experiment_2,
                            _series_id_2,
                            _experiment_2_properties['cells_ids'][_cell_2_id],
                            _cell_2_fiber_densities
                        )

                        _cell_1_fiber_densities_filtered, _cell_2_fiber_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _cell_1_fiber_densities, _cell_2_fiber_densities
                            )

                        # ignore small arrays
                        if len(_cell_1_fiber_densities_filtered) < \
                                MINIMUM_CORRELATION_TIME_FRAMES[_experiment_1]:
                            continue

                        _correlation = compute_lib.correlation(
                            compute_lib.derivative(_cell_1_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_cell_2_fiber_densities_filtered, _n=DERIVATIVE)
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
                'color': '#ea8500'
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
        _filename='plot_real_' + str(_real_cells) + '_offset_z_' + str(_offset_y) + '_high_time_' +
                  str(_high_time_resolution) + '_band_' + str(_band)
    )

    # correct match probability plot
    _y = [_mean_correct_match_probability] * (MAX_RANK - 1) + [_mean_correct_match_probability *
                                                               (_mean_cells_potential_matches - MAX_RANK + 1)]
    _fig = go.Figure(
        data=go.Bar(
            x=_x,
            y=_y,
            marker={
                'color': '#ea8500'
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
        _filename='plot_real_' + str(_real_cells) + '_offset_z_' + str(_offset_y) + '_high_time_' +
                  str(_high_time_resolution) + '_band_' + str(_band) + '_correct_match_prob'
    )


if __name__ == '__main__':
    main()
