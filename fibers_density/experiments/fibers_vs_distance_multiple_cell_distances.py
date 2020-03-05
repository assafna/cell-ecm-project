import os

import numpy as np
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import scatter, save

EXPERIMENTS = ['SN16']
TIME_POINT = 18
CELLS_DISTANCES = [5, 6, 7, 9]
OFFSET_X_STEP = 0.1
BAND = True
OFFSET_Z = 0
OFFSET_Y = 0
OUT_OF_BOUNDARIES = False
OFFSET_X_END = {
    5: 5,
    6: 6,
    7: 7,
    9: 9
}


def main():
    _x_array = []
    _y_array = []
    _names_array = []
    for _distance in CELLS_DISTANCES:
        print('Cells Distance ' + str(_distance))
        _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
        _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
        _experiments = filtering.by_real_cells(_experiments)
        _experiments = filtering.by_distance(_experiments, _distance, _time_point=TIME_POINT - 1)
        if BAND:
            _experiments = filtering.by_band(_experiments)

        _offsets_x = np.arange(start=0, stop=OFFSET_X_END[_distance] + OFFSET_X_STEP, step=OFFSET_X_STEP)

        _arguments = []
        for _tuple in _experiments:
            _experiment, _series_id, _group = _tuple
            for _offset_x in _offsets_x:
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
                    'cell_id': 'left_cell',
                    'direction': 'inside',
                    'time_point': TIME_POINT - 1
                })

        _rois_dictionary, _rois_to_compute = \
            compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x'])
        _fibers_densities = compute.fibers_densities(_rois_to_compute)

        _cells_distance_fibers_densities = [[] for _i in range(len(_offsets_x))]
        for _tuple in _experiments:
            _experiment, _series_id, _group = _tuple
            _offset_index = 0
            _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

            for _offset_x in _offsets_x:
                _fibers_density = _fibers_densities[
                    _rois_dictionary[(_experiment, _series_id, _group, _offset_x)][TIME_POINT - 1]
                ]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _normalized_fibers_density = compute_lib.z_score(
                    _x=_fibers_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                _cells_distance_fibers_densities[_offset_index].append(_normalized_fibers_density)
                _offset_index += 1

        _x_array.append(_offsets_x)
        _y_array.append(_cells_distance_fibers_densities)
        _names_array.append('Cells Distance ' + str(_distance))

    _fig = scatter.create_error_bars_plot(
        _x_array=_x_array,
        _y_array=_y_array,
        _names_array=_names_array,
        _modes_array=['lines+markers'] * len(CELLS_DISTANCES),
        _dashes_array=['solid'] * len(CELLS_DISTANCES),
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
