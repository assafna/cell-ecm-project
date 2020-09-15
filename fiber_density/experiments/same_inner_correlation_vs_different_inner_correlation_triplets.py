import os
import sys
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

MINIMUM_TIME_FRAMES = 18
MINIMUM_PAIR_DISTANCE = 4

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0


def main():
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    _tuples = filtering.by_time_frames_amount(_tuples, MINIMUM_TIME_FRAMES)
    _tuples = filtering.by_pair_distance_range(_tuples, [MINIMUM_PAIR_DISTANCE, sys.maxsize])
    _triplets = filtering.by_triplets(_tuples)
    print('Total triplets:', len(_triplets))

    _arguments = []
    for _triplet in _triplets:
        for _tuple in _triplet:
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
                    'offset_y': OFFSET_Y,
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

    _same_correlations_arrays = [[] for _i in _triplets]
    _different_correlations_arrays = [[] for _i in _triplets]
    _names_array = []
    for _triplet_index, _triplet in enumerate(_triplets):
        for _same_index in tqdm(range(len(_triplet)), desc='Main loop'):
            _same_tuple = _triplet[_same_index]
            _same_experiment, _same_series, _same_group = _same_tuple

            _same_left_cell_fiber_densities = \
                _experiments_fiber_densities[
                    (_same_experiment, _same_series, _same_group, 'left_cell')
                ]
            _same_right_cell_fiber_densities = \
                _experiments_fiber_densities[
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

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_triplet)):
                if _same_index != _different_index:
                    _different_tuple = _triplet[_different_index]
                    _different_experiment, _different_series, _different_group = \
                        _different_tuple
                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fiber_densities = _experiments_fiber_densities[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _same_cell_id
                        )]
                        _different_fiber_densities = _experiments_fiber_densities[(
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

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fiber_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fiber_densities_filtered, _n=DERIVATIVE)
                        )

                        _same_correlations_arrays[_triplet_index].append(_same_correlation)
                        _different_correlations_arrays[_triplet_index].append(_different_correlation)

        _names_array.append('Triplet #' + str(_triplet_index + 1))

    print('Total points:', len(np.array(_same_correlations_arrays).flatten()))
    _same_minus_different = \
        np.array(_same_correlations_arrays).flatten() - np.array(_different_correlations_arrays).flatten()
    print('Wilcoxon of same minus different around the zero:')
    print(wilcoxon(_same_minus_different))
    print('Higher same amount:', (_same_minus_different > 0).sum() /
          len(_same_minus_different))

    # plot
    _colors_array = ['green', 'blue', config.colors(1)]
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_same_correlations_array,
                y=_different_correlations_array,
                name=_name,
                mode='markers',
                marker={
                    'size': 15,
                    'color': _color
                },
                opacity=0.7
            ) for _same_correlations_array, _different_correlations_array, _name, _color in
            zip(_same_correlations_arrays, _different_correlations_arrays, _names_array, _colors_array)
        ],
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
            'legend': {
                'xanchor': 'left',
                'x': 0.1,
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white'
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
