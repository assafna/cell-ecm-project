import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import pearsonr

from libs import compute_lib
from libs.experiments import load, compute, filtering, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

EXPERIMENTS = ['SN16']
BAND = True
TIME_POINT = 18
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_real_cells(_experiments)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    print(len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _cells_distance = \
            compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point=TIME_POINT - 1)
        if _cells_distance / 2 - 0.5 - ROI_LENGTH >= OFFSET_X:
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
                    'time_point': TIME_POINT - 1
                })
        else:
            print(_tuple)

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _x_array = []
    _y_array = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _cell_id in ['left_cell', 'right_cell']:
            if (_experiment, _series_id, _group, _cell_id) in _rois_dictionary:
                _cells_distance = compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group,
                                                                                 _time_point=0)
                _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _cell_id)][0]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    print('Out of boundaries', _tuple)
                    continue

                _normalized_fibers_density = compute_lib.z_score(
                    _x=_fibers_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                if not np.isnan(_normalized_fibers_density):
                    _x_array.append(_cells_distance)
                    _y_array.append(_normalized_fibers_density)
                else:
                    print('nan', _tuple)

    print(len(_x_array))
    print(pearsonr(_x_array, _y_array))

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
                'title': 'Cells Distance'
            },
            'yaxis': {
                'title': 'Fibers Density Z-score'
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_offset_x_' + str(OFFSET_X)
    )


if __name__ == '__main__':
    main()
