import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
ALIGNMENT_OFFSET_Y = 0
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]

START_TIME_FRAME_Z_SCORE = 5
BLOCK_ALREADY_DENSE = False


def align_by_z_score(_experiments, _experiments_fiber_densities):
    _experiments_fiber_densities_aligned = {}
    for _tuple in tqdm(_experiments, desc='Temporal aligning experiments'):
        _experiment, _series_id, _group = _tuple

        _left_cell_alignment_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', ALIGNMENT_OFFSET_Y)]
        _right_cell_alignment_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell', ALIGNMENT_OFFSET_Y)]

        _left_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', OFFSET_Y)]
        _right_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell', OFFSET_Y)]

        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        _left_cell_value_aligned, _right_cell_value_aligned = None, None
        for _time_frame, (_left_cell_fiber_density, _right_cell_fiber_density) in \
                enumerate(zip(_left_cell_alignment_fiber_densities, _right_cell_alignment_fiber_densities)):
            _normalized_left_cell_fiber_density = compute_lib.z_score(
                _x=_left_cell_fiber_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )
            _normalized_right_cell_fiber_density = compute_lib.z_score(
                _x=_right_cell_fiber_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )
            _mean_z_score = (_normalized_left_cell_fiber_density + _normalized_right_cell_fiber_density) / 2
            if _mean_z_score > START_TIME_FRAME_Z_SCORE:
                if _time_frame == 0 and BLOCK_ALREADY_DENSE:
                    _left_cell_value_aligned = []
                    _right_cell_value_aligned = []
                else:
                    _left_cell_value_aligned = _left_cell_fiber_densities[_time_frame:]
                    _right_cell_value_aligned = _right_cell_fiber_densities[_time_frame:]
                break

        # in case z-score not reached
        if _left_cell_value_aligned is None:
            _left_cell_value_aligned = []
            _right_cell_value_aligned = []

        _experiments_fiber_densities_aligned[(_experiment, _series_id, _group, 'left_cell')] = \
            _left_cell_value_aligned
        _experiments_fiber_densities_aligned[(_experiment, _series_id, _group, 'right_cell')] = \
            _right_cell_value_aligned

    return _experiments_fiber_densities_aligned


def main(_high_temporal_resolution=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        for _cell_id in ['left_cell', 'right_cell']:
            _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })
            if ALIGNMENT_OFFSET_Y != OFFSET_Y:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                    'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                    'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                    'offset_x': OFFSET_X,
                    'offset_y': ALIGNMENT_OFFSET_Y,
                    'offset_z': OFFSET_Z,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_points': _latest_time_frame
                })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'offset_y'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _experiments_fiber_densities_aligned = align_by_z_score(_tuples, _experiments_fiber_densities)
    _tuples_by_experiment = organize.by_experiment(_tuples)

    _same_correlations_array = []
    _different_correlations_array = []
    _valid_tuples = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _same_index in tqdm(range(len(_experiment_tuples)), desc='Main loop'):
            _same_tuple = _experiment_tuples[_same_index]
            _same_experiment, _same_series, _same_group = _same_tuple

            _same_left_cell_fiber_densities = \
                _experiments_fiber_densities_aligned[
                    (_same_experiment, _same_series, _same_group, 'left_cell')
                ]
            _same_right_cell_fiber_densities = \
                _experiments_fiber_densities_aligned[
                    (_same_experiment, _same_series, _same_group, 'right_cell')
                ]

            _same_properties = \
                load.group_properties(_same_experiment, _same_series, _same_group)
            _same_left_cell_fiber_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['left_cell'],
                _same_left_cell_fiber_densities
            )
            _same_right_cell_fiber_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['right_cell'],
                _same_right_cell_fiber_densities
            )

            _same_left_cell_fiber_densities_filtered, _same_right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fiber_densities, _same_right_cell_fiber_densities
                )

            # ignore small arrays
            if len(_same_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_same_experiment):
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_experiment_tuples)):
                if _same_index != _different_index:
                    _different_tuple = _experiment_tuples[_different_index]
                    _different_experiment, _different_series, _different_group = \
                        _different_tuple
                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fiber_densities = _experiments_fiber_densities_aligned[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _same_cell_id
                        )]
                        _different_fiber_densities = _experiments_fiber_densities_aligned[(
                            _different_experiment,
                            _different_series,
                            _different_group,
                            _different_cell_id
                        )]

                        _different_properties = load.group_properties(
                            _different_experiment, _different_series, _different_group
                        )
                        _same_fiber_densities = compute.remove_blacklist(
                            _same_experiment,
                            _same_series,
                            _same_properties['cells_ids'][_same_cell_id],
                            _same_fiber_densities
                        )
                        _different_fiber_densities = compute.remove_blacklist(
                            _different_experiment,
                            _different_series,
                            _different_properties['cells_ids'][_different_cell_id],
                            _different_fiber_densities
                        )

                        _same_fiber_densities_filtered, _different_fiber_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fiber_densities, _different_fiber_densities
                            )

                        # ignore small arrays
                        if len(_same_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_different_experiment):
                            continue

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                        )

                        _same_correlations_array.append(_same_correlation)
                        _different_correlations_array.append(_different_correlation)

                        if _same_tuple not in _valid_tuples:
                            _valid_tuples.append(_same_tuple)

    print('Total tuples:', len(_valid_tuples))
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
                'color': '#ea8500'
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
        _filename='plot_high_time_' + str(_high_temporal_resolution)
    )


if __name__ == '__main__':
    main()
