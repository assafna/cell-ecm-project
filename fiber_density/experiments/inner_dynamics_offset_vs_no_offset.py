import os
from itertools import product

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import filtering, load, compute, paths
from libs.experiments.config import all_experiments, QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

OFFSET_X = 0
OFFSET_Y = 1
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]


def main(_high_temporal_resolution=False, _pair_distance_range=None):
    if _pair_distance_range is None:
        _pair_distance_range = PAIR_DISTANCE_RANGE

    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )
    _experiments = ['SN16']

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, _pair_distance_range)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
        for _cell_id, _offset_y in product(['left_cell', 'right_cell'], [0, OFFSET_Y]):
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
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id', 'offset_y'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _no_offset_derivatives_array = []
    _offset_derivatives_array = []
    for _tuple in tqdm(_tuples, desc='Main loop'):
        _experiment, _series_id, _group = _tuple

        _no_offset_left_cell_fiber_densities = _experiments_fiber_densities[
            (_experiment, _series_id, _group, 'left_cell', 0)]
        _no_offset_right_cell_fiber_densities = _experiments_fiber_densities[
            (_experiment, _series_id, _group, 'right_cell', 0)]
        _offset_left_cell_fiber_densities = _experiments_fiber_densities[
            (_experiment, _series_id, _group, 'left_cell', OFFSET_Y)]
        _offset_right_cell_fiber_densities = _experiments_fiber_densities[
            (_experiment, _series_id, _group, 'right_cell', OFFSET_Y)]

        _properties = \
            load.group_properties(_experiment, _series_id, _group)
        _no_offset_left_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series_id,
            _properties['cells_ids']['left_cell'],
            _no_offset_left_cell_fiber_densities
        )
        _no_offset_right_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series_id,
            _properties['cells_ids']['right_cell'],
            _no_offset_right_cell_fiber_densities
        )
        _offset_left_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series_id,
            _properties['cells_ids']['left_cell'],
            _offset_left_cell_fiber_densities
        )
        _offset_right_cell_fiber_densities = compute.remove_blacklist(
            _experiment,
            _series_id,
            _properties['cells_ids']['right_cell'],
            _offset_right_cell_fiber_densities
        )

        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        for _cell_pair in [(_no_offset_left_cell_fiber_densities, _offset_left_cell_fiber_densities),
                           (_no_offset_right_cell_fiber_densities, _offset_right_cell_fiber_densities)]:
            for _time_frame in range(1, min(len(_cell_pair[0]), len(_cell_pair[1]))):
                _no_offset_cell_fiber_density_previous = _cell_pair[0][_time_frame - 1]
                _no_offset_cell_fiber_density = _cell_pair[0][_time_frame]
                _offset_cell_fiber_density_previous = _cell_pair[1][_time_frame - 1]
                _offset_cell_fiber_density = _cell_pair[1][_time_frame]

                if any([_no_offset_cell_fiber_density_previous[1], _no_offset_cell_fiber_density[1],
                        _offset_cell_fiber_density_previous[1], _offset_cell_fiber_density[1]]):
                    continue

                _no_offset_cell_fiber_density_previous_z_score = compute_lib.z_score(
                    _x=_no_offset_cell_fiber_density_previous[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )
                _no_offset_cell_fiber_density_z_score = compute_lib.z_score(
                    _x=_no_offset_cell_fiber_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )
                _offset_cell_fiber_density_previous_z_score = compute_lib.z_score(
                    _x=_offset_cell_fiber_density_previous[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )
                _offset_cell_fiber_density_z_score = compute_lib.z_score(
                    _x=_offset_cell_fiber_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                _no_offset_derivative = _no_offset_cell_fiber_density_z_score - _no_offset_cell_fiber_density_previous_z_score
                _offset_derivative = _offset_cell_fiber_density_z_score - _offset_cell_fiber_density_previous_z_score

                _no_offset_derivatives_array.append(_no_offset_derivative)
                _offset_derivatives_array.append(_offset_derivative)

    print('Total points:', len(_no_offset_derivatives_array))
    _offset_minus_no_offset = \
        np.array(_offset_derivatives_array) - np.array(_no_offset_derivatives_array)
    print('Wilcoxon of offset minus no offset around the zero:')
    print(wilcoxon(_offset_minus_no_offset))
    print('Higher offset amount:', (_offset_minus_no_offset > 0).sum() /
          len(_offset_minus_no_offset))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_no_offset_derivatives_array,
            y=_offset_derivatives_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'No offset derivative',
                # 'zeroline': False,
                # 'range': [-1.1, 1.2],
                # 'tickmode': 'array',
                # 'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Offset derivative',
                # 'zeroline': False,
                # 'range': [-1.1, 1.2],
                # 'tickmode': 'array',
                # 'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'shapes': [
                # {
                #     'type': 'line',
                #     'x0': -1,
                #     'y0': -1,
                #     'x1': -1,
                #     'y1': 1,
                #     'line': {
                #         'color': 'black',
                #         'width': 2
                #     }
                # },
                # {
                #     'type': 'line',
                #     'x0': -1,
                #     'y0': -1,
                #     'x1': 1,
                #     'y1': -1,
                #     'line': {
                #         'color': 'black',
                #         'width': 2
                #     }
                # },
                {
                    'type': 'line',
                    'x0': -5,
                    'y0': -5,
                    'x1': 5,
                    'y1': 5,
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
