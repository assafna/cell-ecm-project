import os
from itertools import product

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, compute, filtering, organize, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

PAIR_DISTANCE = 7

# according to pair distance
OFFSET_X_END = {
    5: 2.8,
    7: 5.3
}

OFFSET_X_STEP = 0.1
OFFSETS_X = np.arange(start=0, stop=OFFSET_X_END[PAIR_DISTANCE] + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0
OFFSET_Z = 0


def main():
    print('Single Cell')
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=True,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.minimum_time_frames_for_correlation(_experiments[0]))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.minimum_time_frames_for_correlation(_experiment)
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right']):
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': _time_frame - 1
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'direction'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _tuples = organize.by_single_cell_id(_tuples)

    _single_cell_fiber_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in _tuples:
        _experiment, _series_id, _cell_id = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Cell ID:', _cell_id, sep='\t')
        _offset_index = 0
        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        for _offset_x in OFFSETS_X:
            _cell_fiber_densities = []
            for _cell_tuple in _tuples[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _offset_x, _direction)][0]
                    _fiber_density = _fiber_densities[_window_tuple]

                    if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                        continue

                    _normalized_fiber_density = compute_lib.z_score(
                        _x=_fiber_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    _cell_fiber_densities.append(_normalized_fiber_density)

            if len(_cell_fiber_densities) > 0:
                _single_cell_fiber_densities[_offset_index].append(np.mean(_cell_fiber_densities))
            _offset_index += 1

    print('Pairs')
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.minimum_time_frames_for_correlation(_experiments[0]))
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_pair_distance(_tuples, PAIR_DISTANCE)
    _tuples = filtering.by_band(_tuples)

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.minimum_time_frames_for_correlation(_experiment)
        for _offset_x in OFFSETS_X:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'left_cell',
                'direction': 'inside',
                'time_point': _time_frame - 1
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _pairs_fiber_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Group:', _group, sep='\t')
        _offset_index = 0
        _normalization = load.normalization_series_file_data(_experiment, _series_id)

        # take offsets based on pair distance
        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_coordinates = [list(_properties['time_points'][0]['left_cell']['coordinates'].values())]
        _right_cell_coordinates = [list(_properties['time_points'][0]['right_cell']['coordinates'].values())]
        _pair_distance = compute.pair_distance_in_cell_size(
            _experiment, _series_id, _left_cell_coordinates, _right_cell_coordinates)
        _edges_distance = _pair_distance - 1
        _max_x_offset = _edges_distance - QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER

        for _offset_x in OFFSETS_X:
            if _offset_x > _max_x_offset:
                break

            _fiber_density = _fiber_densities[_windows_dictionary[(_experiment, _series_id, _group, _offset_x)][0]]

            if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                continue

            _normalized_fiber_density = compute_lib.z_score(
                _x=_fiber_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )

            _pairs_fiber_densities[_offset_index].append(_normalized_fiber_density)
            _offset_index += 1

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _pairs_fiber_densities],
                name='Pairs',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _pairs_fiber_densities],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'solid'}
            ),
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _single_cell_fiber_densities],
                name='Single Cell',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _single_cell_fiber_densities],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'dash'}
            )
        ],
        layout={
            'xaxis_title': 'Distance from left cell (cell size)',
            'yaxis_title': 'Fiber density (z-score)'
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_distance_' + str(PAIR_DISTANCE)
    )


if __name__ == '__main__':
    main()
