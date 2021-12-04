import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import filtering, load, compute, paths
from libs.experiments.config import all_experiments, QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]


def main(_high_temporal_resolution=True, _offset_y=0.5):
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
    _tuples = filtering.by_real_pairs(_tuples)
    # _tuples = filtering.by_band(_tuples)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
        for _cell_id, _direction in product(['left_cell', 'right_cell'], ['inside', 'outside']):
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

    _inner_correlations = []
    _outer_correlations = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _left_inner_tuple = (_experiment, _series_id, _group, 'left_cell', 'inside')
        _left_outer_tuple = (_experiment, _series_id, _group, 'left_cell', 'outside')
        _right_inner_tuple = (_experiment, _series_id, _group, 'right_cell', 'inside')
        _right_outer_tuple = (_experiment, _series_id, _group, 'right_cell', 'outside')

        if _left_inner_tuple not in _windows_dictionary or _left_outer_tuple not in _windows_dictionary or \
                _right_inner_tuple not in _windows_dictionary or _right_outer_tuple not in _windows_dictionary:
            continue

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_inner_fiber_densities = _experiments_fiber_densities[_left_inner_tuple]
        _left_outer_fiber_densities = _experiments_fiber_densities[_left_outer_tuple]
        _right_inner_fiber_densities = _experiments_fiber_densities[_right_inner_tuple]
        _right_outer_fiber_densities = _experiments_fiber_densities[_right_outer_tuple]

        _left_inner_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_inner_fiber_densities)
        _left_outer_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_outer_fiber_densities)
        _right_inner_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_inner_fiber_densities)
        _right_outer_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_outer_fiber_densities)

        _left_inner_fiber_densities, _right_inner_fiber_densities = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_inner_fiber_densities, _right_inner_fiber_densities)
        _left_outer_fiber_densities, _right_outer_fiber_densities = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_outer_fiber_densities, _right_outer_fiber_densities)

        # ignore small arrays
        _minimum_time_frame_for_correlation = compute.minimum_time_frames_for_correlation(_experiment)
        if len(_left_inner_fiber_densities) < _minimum_time_frame_for_correlation or \
                len(_left_outer_fiber_densities) < _minimum_time_frame_for_correlation:
            continue

        _inner_correlations.append(
            compute_lib.correlation(
                compute_lib.derivative(_left_inner_fiber_densities, _n=DERIVATIVE),
                compute_lib.derivative(_right_inner_fiber_densities, _n=DERIVATIVE)))
        _outer_correlations.append(
            compute_lib.correlation(
                compute_lib.derivative(_left_outer_fiber_densities, _n=DERIVATIVE),
                compute_lib.derivative(_right_outer_fiber_densities, _n=DERIVATIVE)))

    print('Total pairs:', len(_inner_correlations))
    print('Pearson correlation:', compute_lib.correlation(_inner_correlations, _outer_correlations, _with_p_value=True))
    print('Wilcoxon of inner around the zero:', wilcoxon(_inner_correlations))
    print('Wilcoxon of outer around the zero:', wilcoxon(_outer_correlations))
    _inner_minus_outer = np.array(_inner_correlations) - np.array(_outer_correlations)
    print('Wilcoxon of inner minus outer:', wilcoxon(_inner_minus_outer))
    print('Higher inner amount:', (_inner_minus_outer > 0).sum() / len(_inner_minus_outer))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_inner_correlations,
            y=_outer_correlations,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Inner correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Outer correlation',
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
        _filename='plot_high_' + str(_high_temporal_resolution) + '_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
