import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

EXPERIMENTS = ['SN16']
TIME_POINT = 18
# CELLS_DISTANCES = [7]
OFFSET_X_STEP = 0.2
BAND = True
OFFSET_Z = 0
OFFSET_Y = 0
OUT_OF_BOUNDARIES = False


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_real_cells(_experiments)
    # _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    _max_offsets_x = []
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _cells_distance = \
            compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point=TIME_POINT - 1)
        _offsets_x = \
            np.arange(start=0, stop=_cells_distance / 2 - 0.5 - ROI_LENGTH + OFFSET_X_STEP, step=OFFSET_X_STEP)
        if len(_offsets_x) > len(_max_offsets_x):
            _max_offsets_x = _offsets_x
        for _offset_x in _offsets_x:
            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': ROI_LENGTH,
                    'length_y': ROI_HEIGHT,
                    'length_z': ROI_WIDTH,
                    'offset_x': _offset_x,
                    'offset_y': OFFSET_Y,
                    'offset_z': OFFSET_Z,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': TIME_POINT - 1
                })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _x_array = []
    _y_array = []
    for _offset_x in _max_offsets_x:
        _cells_distances = []
        _z_scores = []
        for _tuple in _experiments:
            _experiment, _series_id, _group = _tuple
            for _cell_id in ['left_cell', 'right_cell']:
                if (_experiment, _series_id, _group, _offset_x, _cell_id) in _rois_dictionary:
                    _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
                    _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _offset_x, _cell_id)][0]
                    _fibers_density = _fibers_densities[_roi_tuple]

                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    if not np.isnan(_normalized_fibers_density):
                        _cells_distances.append(
                            compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group,
                                                                           _time_point=0)
                        )
                        _z_scores.append(_normalized_fibers_density)
        _x_array.append(_offset_x)
        _y_array.append(compute_lib.correlation(_cells_distances, _z_scores))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_x_array,
            y=_y_array,
            mode='markers',
            marker={
                'size': 15
            },
            opacity=0.7
        ),
        layout={
            'xaxis': {
                'title': 'Distance from Cell (cell size)',
                # 'zeroline': False
            },
            'yaxis': {
                'title': 'Correlation',
                # 'range': [-1.7, 16],
                # 'zeroline': False,
                # 'tickmode': 'array',
                # 'tickvals': [0, 4, 8, 12]
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
