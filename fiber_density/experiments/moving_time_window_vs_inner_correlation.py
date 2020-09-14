import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41', 'SN44', 'SN45']
}
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 1
TIME_RESOLUTION = {
    False: 15,
    True: 5
}
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
BAND = True
DIRECTION = 'inside'
MOVING_WINDOW_LENGTH = {
    False: 10,
    True: 50
}
TIME_FRAME_STEP = {
    False: 1,
    True: 10
}
END_TIME_FRAME = {
    False: 20,
    True: 200
}


def main(_high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = \
        filtering.by_time_frames_amount(_experiments, _time_frames=MOVING_WINDOW_LENGTH[_high_time_resolution])
    _experiments = filtering.by_pair_distance_range(_experiments, _distance_range=PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_real_pairs(_experiments, _real_pairs=REAL_CELLS)
    _experiments = filtering.by_fake_static_pairs(_experiments, _fake_static_pairs=STATIC)
    _experiments = filtering.by_band(_experiments, _band=BAND)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_frame = len(_properties['time_points'])
        if DIRECTION == 'inside':
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
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    for _tuple in _experiments:
        _correlations = []
        _experiment, _series_id, _group = _tuple

        _left_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        for _start_time_frame in \
                range(0, END_TIME_FRAME[_high_time_resolution], TIME_FRAME_STEP[_high_time_resolution]):

            _left_cell_fiber_densities_window = _left_cell_fiber_densities[
                                                _start_time_frame:_start_time_frame +
                                                MOVING_WINDOW_LENGTH[_high_time_resolution]]
            _right_cell_fiber_densities_window = _right_cell_fiber_densities[
                                                 _start_time_frame:_start_time_frame +
                                                 MOVING_WINDOW_LENGTH[_high_time_resolution]]

            _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_fiber_densities_window, _right_cell_fiber_densities_window)

            # ignore small arrays
            if len(_left_cell_fiber_densities_filtered) < MOVING_WINDOW_LENGTH[_high_time_resolution]:
                _correlations.append(None)
                continue

            _correlations.append(compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            ))

        # plot
        _fig = go.Figure(
            data=go.Scatter(
                x=np.arange(
                    start=0,
                    stop=len(_correlations),
                    step=1) * TIME_RESOLUTION[_high_time_resolution] * TIME_FRAME_STEP[_high_time_resolution],
                y=_correlations,
                mode='lines+markers',
                line={'dash': 'solid'}
            ),
            layout={
                'xaxis': {
                    'title': 'Window start time (minutes)',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Inner correlation',
                    'zeroline': False
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_' + str(_experiment) + '_' + str(_series_id) + '_' + str(_group)
        )


if __name__ == '__main__':
    main()
