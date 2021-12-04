import os
from itertools import product

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def main(_directions=None):
    if _directions is None:
        _directions = ['inside', 'outside']

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
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
        for _cell_id, _direction in product(['left_cell', 'right_cell'], _directions):
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
                'direction': _direction,
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'direction'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    for _direction in _directions:
        _y_arrays = [[] for _i in DERIVATIVES]
        for _tuple in tqdm(_tuples, desc='Experiments loop'):
            _experiment, _series_id, _group = _tuple

            if (_experiment, _series_id, _group, 'left_cell', _direction) not in _windows_dictionary or \
                    (_experiment, _series_id, _group, 'right_cell', _direction) not in _windows_dictionary:
                continue

            _properties = load.group_properties(_experiment, _series_id, _group)

            _left_cell_fiber_densities = \
                _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', _direction)]
            _right_cell_fiber_densities = \
                _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell', _direction)]

            _left_cell_fiber_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
            _right_cell_fiber_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

            if not OUT_OF_BOUNDARIES:
                _left_cell_fiber_densities, _right_cell_fiber_densities = \
                    compute.longest_same_indices_shared_in_borders_sub_array(
                        _left_cell_fiber_densities, _right_cell_fiber_densities
                    )
            else:
                _left_cell_fiber_densities = [_fiber_density[0] for _fiber_density in _left_cell_fiber_densities]
                _right_cell_fiber_densities = [_fiber_density[0] for _fiber_density in _right_cell_fiber_densities]

            # ignore small arrays
            _minimum_time_frame_for_correlation = compute.minimum_time_frames_for_correlation(_experiment)
            if len(_left_cell_fiber_densities) < _minimum_time_frame_for_correlation or \
                    len(_right_cell_fiber_densities) < _minimum_time_frame_for_correlation:
                continue

            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fiber_densities, _n=_derivative),
                    compute_lib.derivative(_right_cell_fiber_densities, _n=_derivative)
                ))

        print('Direction:', _direction)
        print('Total pairs:', len(_y_arrays[0]))
        print('Wilcoxon around the zero')
        for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
            print('Derivative:', _derivative, wilcoxon(_y_array))

        # plot
        _y_title = 'Inner correlation' if _direction == 'inside' else 'Outer correlation'
        _colors_array = config.colors(3)
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_derivative,
                    boxpoints='all',
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
                    opacity=0.7,
                    showlegend=False
                ) for _y, _derivative, _color in zip(_y_arrays, DERIVATIVES_TEXT, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Fiber density derivative',
                    'zeroline': False
                },
                'yaxis': {
                    'title': _y_title,
                    'range': [-1, 1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_direction_' + _direction
        )


if __name__ == '__main__':
    main()
