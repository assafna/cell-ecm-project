import os
import sys
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

EXPERIMENTS = ['SN16']
MINIMUM_TIME_POINTS = 18
MINIMUM_CELLS_DISTANCE = 4
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_real_cells(_experiments)
    _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance_range(_experiments, [MINIMUM_CELLS_DISTANCE, sys.maxsize])
    _triplets = filtering.by_triplets(_experiments)
    print('Total triplets:', len(_triplets))

    _arguments = []
    for _triplet in _triplets:
        for _tuple in _triplet:
            _experiment, _series_id, _group = _tuple

            # stop when windows are overlapping
            _properties = load.group_properties(_experiment, _series_id, _group)
            _latest_time_point = len(_properties['time_points'])
            for _time_point in range(len(_properties['time_points'])):
                _cells_distance = \
                    compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
                if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                    _latest_time_point = _time_point - 1
                    break

            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': ROI_LENGTH,
                    'length_y': ROI_HEIGHT,
                    'length_z': ROI_WIDTH,
                    'offset_x': OFFSET_X,
                    'offset_y': OFFSET_Y,
                    'offset_z': OFFSET_Z,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_points': _latest_time_point
                })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _same_correlations_arrays = [[] for _i in _triplets]
    _different_correlations_arrays = [[] for _i in _triplets]
    _names_array = []
    for _triplet_index, _triplet in enumerate(_triplets):
        for _same_index in tqdm(range(len(_triplet)), desc='Main loop'):
            _same_tuple = _triplet[_same_index]
            _same_experiment, _same_series, _same_group = _same_tuple

            _same_left_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_same_experiment, _same_series, _same_group, 'left_cell')
                ]
            _same_right_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_same_experiment, _same_series, _same_group, 'right_cell')
                ]

            _same_properties = \
                load.group_properties(_same_experiment, _same_series, _same_group)
            _same_left_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['left_cell'],
                _same_left_cell_fibers_densities
            )
            _same_right_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['right_cell'],
                _same_right_cell_fibers_densities
            )

            _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fibers_densities, _same_right_cell_fibers_densities
                )

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_triplet)):
                if _same_index != _different_index:
                    _different_tuple = _triplet[_different_index]
                    _different_experiment, _different_series, _different_group = \
                        _different_tuple
                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fibers_densities = _experiments_fibers_densities[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _same_cell_id
                        )]
                        _different_fibers_densities = _experiments_fibers_densities[(
                            _different_experiment,
                            _different_series,
                            _different_group,
                            _different_cell_id
                        )]

                        _different_properties = load.group_properties(
                            _different_experiment, _different_series, _different_group
                        )
                        _same_fibers_densities = compute.remove_blacklist(
                            _same_experiment,
                            _same_series,
                            _same_properties['cells_ids'][_same_cell_id],
                            _same_fibers_densities
                        )
                        _different_fibers_densities = compute.remove_blacklist(
                            _different_experiment,
                            _different_series,
                            _different_properties['cells_ids'][_different_cell_id],
                            _different_fibers_densities
                        )

                        _same_fibers_densities_filtered, _different_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fibers_densities, _different_fibers_densities
                            )

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
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
    _colors_array = ['green', 'blue', 'red']
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
