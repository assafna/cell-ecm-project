import os

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
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
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
BAND = True
DIRECTION = 'inside'
MOVING_WINDOW_LENGTH = {
    False: 10,
    True: 50
}
TIME_POINT_STEP = {
    False: 1,
    True: 10
}
END_TIME_POINT = {
    False: 20,
    True: 200
}


def main(_high_time_resolution=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = \
        filtering.by_time_points_amount(_experiments, _time_points=MOVING_WINDOW_LENGTH[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, _distance_range=CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=BAND)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
        if DIRECTION == 'inside':
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
                'direction': DIRECTION,
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _correlations_by_start_time_point = []
    for _start_time_point in tqdm(range(0, END_TIME_POINT[_high_time_resolution],
                                        TIME_POINT_STEP[_high_time_resolution]), desc='Time loop'):
        _correlations = []
        for _tuple in _experiments:
            _experiment, _series_id, _group = _tuple

            _left_cell_fibers_densities = \
                _experiments_fibers_densities[(_experiment, _series_id, _group, 'left_cell')]
            _right_cell_fibers_densities = \
                _experiments_fibers_densities[(_experiment, _series_id, _group, 'right_cell')]

            _properties = load.group_properties(_experiment, _series_id, _group)
            _left_cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fibers_densities)
            _right_cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fibers_densities)

            _left_cell_fibers_densities_window = _left_cell_fibers_densities[
                                                 _start_time_point:_start_time_point +
                                                 MOVING_WINDOW_LENGTH[_high_time_resolution]]
            _right_cell_fibers_densities_window = _right_cell_fibers_densities[
                                                  _start_time_point:_start_time_point +
                                                  MOVING_WINDOW_LENGTH[_high_time_resolution]]

            _left_cell_fibers_densities_filtered, _right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_fibers_densities_window, _right_cell_fibers_densities_window)

            # ignore small arrays
            if len(_left_cell_fibers_densities_filtered) < MOVING_WINDOW_LENGTH[_high_time_resolution]:
                continue

            _correlations.append(compute_lib.correlation(
                compute_lib.derivative(_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            ))

        _correlations_by_start_time_point.append(_correlations)

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=np.arange(
                start=0,
                stop=len(_correlations_by_start_time_point),
                step=1) * TIME_RESOLUTION[_high_time_resolution] * TIME_POINT_STEP[_high_time_resolution],
            y=[np.mean(_array) for _array in _correlations_by_start_time_point],
            error_y={
                'type': 'data',
                'array': [np.std(_array) for _array in _correlations_by_start_time_point],
                'thickness': 1
            },
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
