import os
import random

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]
STATIC = False

MAX_RANK = 7
POTENTIAL_MATCHES = 50


def main(_real_cells=True, _offset_y=0.5, _high_temporal_resolution=False, _band=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=_real_cells)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=STATIC)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
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
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _tuples_by_experiment = organize.by_experiment(_tuples)

    _n = 0
    _cells_ranks = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _pivot_tuple in tqdm(_experiment_tuples, desc='Main loop'):
            _pivot_experiment, _pivot_series_id, _pivot_group = _pivot_tuple
            _pivot_experiment_properties = load.group_properties(_pivot_experiment, _pivot_series_id, _pivot_group)

            for _pivot_cell_id, _pivot_cell_correct_match_cell_id in \
                    zip(['left_cell', 'right_cell'], ['right_cell', 'left_cell']):
                _pivot_cell = (_pivot_experiment, _pivot_series_id, _pivot_group, _pivot_cell_id)
                _pivot_cell_correct_match_cell = (_pivot_experiment, _pivot_series_id, _pivot_group, _pivot_cell_correct_match_cell_id)
                _pivot_cell_fiber_densities = _experiments_fiber_densities[_pivot_cell]
                _pivot_cell_fiber_densities = compute.remove_blacklist(
                    _pivot_experiment,
                    _pivot_series_id,
                    _pivot_experiment_properties['cells_ids'][_pivot_cell_id],
                    _pivot_cell_fiber_densities
                )

                _pivot_cell_correlations = []

                # correct match
                _pivot_cell_correct_match_fiber_densities = _experiments_fiber_densities[_pivot_cell_correct_match_cell]
                _pivot_cell_correct_match_fiber_densities = compute.remove_blacklist(
                    _pivot_experiment,
                    _pivot_series_id,
                    _pivot_experiment_properties['cells_ids'][_pivot_cell_correct_match_cell_id],
                    _pivot_cell_correct_match_fiber_densities
                )
                _pivot_cell_fiber_densities_filtered, _pivot_cell_correct_match_fiber_densities_filtered = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _pivot_cell_fiber_densities, _pivot_cell_correct_match_fiber_densities
                    )
                # ignore small arrays
                if len(_pivot_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(
                        _pivot_experiment):
                    continue
                _correlation = compute_lib.correlation(
                    compute_lib.derivative(_pivot_cell_fiber_densities_filtered, _n=DERIVATIVE),
                    compute_lib.derivative(_pivot_cell_correct_match_fiber_densities_filtered, _n=DERIVATIVE)
                )
                _pivot_cell_correlations.append(_correlation)
                _pivot_cell_correct_match_correlation = _correlation

                # create list of potential cells
                _candidate_tuples = []
                for _candidate_tuple in _experiment_tuples:
                    _candidate_experiment, _candidate_series_id, _candidate_group = _candidate_tuple
                    for _candidate_cell_id in ['left_cell', 'right_cell']:
                        _candidate_cell = (_candidate_experiment, _candidate_series_id, _candidate_group, _candidate_cell_id)

                        # skip if same cell or correct match
                        if _candidate_cell == _pivot_cell or _candidate_cell == _pivot_cell_correct_match_cell:
                            continue

                        _candidate_tuples.append(_candidate_cell)

                # compare with each potential candidate, until reached the maximum or nothing to compare with
                while len(_pivot_cell_correlations) < POTENTIAL_MATCHES and len(_candidate_tuples) > 0:

                    # sample randomly
                    _candidate_cell = random.choice(_candidate_tuples)
                    _candidate_experiment, _candidate_series_id, _candidate_group, _candidate_cell_id = _candidate_cell
                    _candidate_experiment_properties = load.group_properties(_candidate_experiment, _candidate_series_id, _candidate_group)

                    _candidate_cell_fiber_densities = _experiments_fiber_densities[_candidate_cell]
                    _candidate_cell_fiber_densities = compute.remove_blacklist(
                        _candidate_experiment,
                        _candidate_series_id,
                        _candidate_experiment_properties['cells_ids'][_candidate_cell_id],
                        _candidate_cell_fiber_densities
                    )

                    _pivot_cell_fiber_densities_filtered, _candidate_cell_fiber_densities_filtered = \
                        compute.longest_same_indices_shared_in_borders_sub_array(
                            _pivot_cell_fiber_densities, _candidate_cell_fiber_densities
                        )

                    # ignore small arrays
                    if len(_pivot_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_pivot_experiment):
                        _candidate_tuples.remove(_candidate_cell)
                        continue

                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_pivot_cell_fiber_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_candidate_cell_fiber_densities_filtered, _n=DERIVATIVE)
                    )

                    _pivot_cell_correlations.append(_correlation)

                # nothing to compare with
                if len(_pivot_cell_correlations) == 1:
                    continue

                # check matchmaking
                _pivot_cell_correct_match_rank = 1
                for _potential_match_correlation in sorted(_pivot_cell_correlations, reverse=True):
                    if _pivot_cell_correct_match_correlation == _potential_match_correlation:
                        break
                    _pivot_cell_correct_match_rank += 1

                _n += 1
                _cells_ranks.append(_pivot_cell_correct_match_rank)

    # results
    _correct_match_probability = 1 / POTENTIAL_MATCHES
    _first_place_correct_matches = sum([1 for _rank in _cells_ranks if _rank == 1])
    _first_place_fraction = _first_place_correct_matches / _n

    print('Matchmaking results:')
    print('Total cells:', _n)
    print('Correct match probability:', round(_correct_match_probability, 2))
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
                  str(_high_temporal_resolution) + '_band_' + str(_band)
    )

    # correct match probability plot
    _y = [_correct_match_probability] * (MAX_RANK - 1) + [_correct_match_probability * (POTENTIAL_MATCHES - MAX_RANK + 1)]
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
                  str(_high_temporal_resolution) + '_band_' + str(_band) + '_correct_match_prob'
    )


if __name__ == '__main__':
    main()
