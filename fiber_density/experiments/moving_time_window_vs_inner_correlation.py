import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]

# according to temporal resolution
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


def main(_high_temporal_resolution=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, _time_frames=MOVING_WINDOW_LENGTH[_high_temporal_resolution])
    _tuples = filtering.by_pair_distance_range(_tuples, _distance_range=PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
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
                'offset_y': OFFSET_Y,
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

    for _tuple in _tuples:
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
                range(0, END_TIME_FRAME[_high_temporal_resolution], TIME_FRAME_STEP[_high_temporal_resolution]):

            _left_cell_fiber_densities_window = _left_cell_fiber_densities[
                                                _start_time_frame:_start_time_frame +
                                                MOVING_WINDOW_LENGTH[_high_temporal_resolution]]
            _right_cell_fiber_densities_window = _right_cell_fiber_densities[
                                                 _start_time_frame:_start_time_frame +
                                                 MOVING_WINDOW_LENGTH[_high_temporal_resolution]]

            _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_fiber_densities_window, _right_cell_fiber_densities_window)

            # ignore small arrays
            if len(_left_cell_fiber_densities_filtered) < MOVING_WINDOW_LENGTH[_high_temporal_resolution]:
                _correlations.append(None)
                continue

            _correlations.append(compute_lib.correlation(
                compute_lib.derivative(_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
            ))

        # plot
        _temporal_resolution = compute.temporal_resolution_in_minutes(_experiment)
        _fig = go.Figure(
            data=go.Scatter(
                x=np.arange(
                    start=0,
                    stop=len(_correlations),
                    step=1) * _temporal_resolution * TIME_FRAME_STEP[_high_temporal_resolution],
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
