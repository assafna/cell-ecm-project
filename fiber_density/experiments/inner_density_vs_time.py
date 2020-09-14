import os

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, compute, paths, filtering
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments
from plotting import save

TEMPORAL_RESOLUTION = 15

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0

DERIVATIVE = 0


def compute_tuples(_tuples):
    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_frame = len(_properties['time_points'])
        _middle_offsets_x = []
        for _time_frame in range(len(_properties['time_points'])):
            _pair_distance = \
                compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame)
            if _time_frame == 0:
                print(_tuple, 'Pair distance:', _pair_distance)
            _middle_offsets_x.append(
                (_time_frame, (_pair_distance - 1) / 2 - QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER / 2))
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
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame,
                'middle_time_frame': -1
            })

        # middle one
        for (_time_frame, _offset_x) in _middle_offsets_x:
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
                'time_point': _time_frame,
                'middle_time_frame': _time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'middle_time_frame'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _tuples_data = []
    for _tuple in tqdm(_tuples, desc='Main loop'):
        _experiment, _series_id, _group = _tuple

        _left_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', -1)]
        _right_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell', -1)]
        _middle_fiber_densities = \
            [_experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', _time_frame)][0]
             for _time_frame in range(0, len(_left_cell_fiber_densities))]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_cell_fiber_densities, _right_cell_fiber_densities)

        if len(_left_cell_fiber_densities_filtered) == 0:
            continue

        _start_time_frame = 0
        for _left in _left_cell_fiber_densities:
            if _left[0] == _left_cell_fiber_densities_filtered[0]:
                break
            _start_time_frame += 1

        _middle_fiber_densities_filtered = \
            [_fiber_density[0] for _fiber_density in _middle_fiber_densities]
        _middle_fiber_densities_filtered = \
            _middle_fiber_densities_filtered[_start_time_frame:
                                             _start_time_frame + len(_left_cell_fiber_densities_filtered)]

        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        _left_cell_fiber_densities_normalized = compute_lib.z_score_array(
            _array=_left_cell_fiber_densities_filtered,
            _average=_normalization['average'],
            _std=_normalization['std']
        )
        _right_cell_fiber_densities_normalized = compute_lib.z_score_array(
            _array=_right_cell_fiber_densities_filtered,
            _average=_normalization['average'],
            _std=_normalization['std']
        )
        _middle_fiber_densities_normalized = compute_lib.z_score_array(
            _array=_middle_fiber_densities_filtered,
            _average=_normalization['average'],
            _std=_normalization['std']
        )

        _left_cell_fiber_densities_normalized = \
            compute_lib.derivative(_left_cell_fiber_densities_normalized, _n=DERIVATIVE)
        _right_cell_fiber_densities_normalized = \
            compute_lib.derivative(_right_cell_fiber_densities_normalized, _n=DERIVATIVE)
        _middle_fiber_densities_normalized = \
            compute_lib.derivative(_middle_fiber_densities_normalized, _n=DERIVATIVE)

        _correlation = \
            compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities_normalized, _n=1),
                compute_lib.derivative(_right_cell_fiber_densities_normalized, _n=1)
            )

        _tuples_data.append([
            _tuple,
            _start_time_frame,
            (
                _left_cell_fiber_densities_normalized,
                _right_cell_fiber_densities_normalized,
                _middle_fiber_densities_normalized
            ),
            _correlation
        ])

    # plots
    _names_array = ['Left cell', 'Right cell', 'Middle']
    _colors_array = ['#844b00', '#ea8500', '#edbc80']
    for _tuple_data in _tuples_data:
        _tuple, _start_time_frame, _y_arrays, _correlation = _tuple_data
        _fig = go.Figure(
            data=[
                go.Scatter(
                    x=np.arange(
                        start=_start_time_frame,
                        stop=_start_time_frame + len(_y_arrays[0]) - DERIVATIVE,
                        step=1) * TEMPORAL_RESOLUTION,
                    y=_y,
                    name=_name,
                    mode='lines',
                    line={
                        'color': _color
                    }
                ) for _y, _name, _color in zip(_y_arrays, _names_array, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Time (minutes)',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Fiber density (z-score)',
                    'zeroline': False,
                    'range': [-3, 8],
                    'tickmode': 'array',
                    'tickvals': [-3, 0, 3, 6]
                },
                'legend': {
                    'xanchor': 'left',
                    'x': 0.1,
                    'yanchor': 'top',
                    'bordercolor': 'black',
                    'borderwidth': 2,
                    'bgcolor': 'white'
                },
            }
        )

        _experiment, _series_id, _group = _tuple
        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_' + _experiment + '_' + str(_series_id) + '_' + _group
        )

        print('Tuple:', _tuple, 'Correlation:', _correlation)


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
    print('Total tuples:', len(_tuples))

    compute_tuples(_tuples)


if __name__ == '__main__':
    main()
