import math
import os
from itertools import product

import numpy as np
from multiprocess.pool import Pool
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import USE_MULTIPROCESSING, CPUS_TO_USE
from libs.experiments import load, compute, config, filtering, organize, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

PAIRS_EXPERIMENTS = ['SN16']
PAIRS_CELLS_DISTANCE = 7
PAIRS_BAND = True
TIME_POINT = 18
OFFSET_X_END = {
    5: 2.8,
    7: 5.3
}
OFFSET_X_STEP = 0.1
OFFSETS_X = np.arange(start=0, stop=OFFSET_X_END[PAIRS_CELLS_DISTANCE] + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False


def main():
    # single cell
    print('Single Cell')
    _experiments = load.experiments_groups_as_tuples(config.SINGLE_CELL)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right']):
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
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': TIME_POINT - 1
            })

    _rois = compute.rois(_arguments)
    _fibers_densities = compute.fibers_densities(_rois)

    _experiments = organize.by_single_cell_id(_experiments)

    _single_cell_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in _experiments:
        _experiment, _series_id, _cell_id = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Cell ID:', _cell_id, sep='\t')
        _offset_index = 0
        _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
        for _offset_x in OFFSETS_X:
            _cell_fibers_densities = []
            for _cell_tuple in _experiments[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _arguments = {
                        'experiment': _experiment,
                        'series_id': _series_id,
                        'group': _group,
                        'length_x': ROI_LENGTH,
                        'length_y': ROI_HEIGHT,
                        'length_z': ROI_WIDTH,
                        'offset_x': _offset_x,
                        'offset_y': OFFSET_Y,
                        'offset_z': OFFSET_Z,
                        'cell_id': 'cell',
                        'direction': _direction,
                        'time_point': TIME_POINT - 1
                    }
                    _roi = compute.roi_time_point(_arguments)
                    _fibers_density = _fibers_densities[_roi]

                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    _cell_fibers_densities.append(_normalized_fibers_density)

            if len(_cell_fibers_densities) > 0:
                _single_cell_fibers_densities[_offset_index].append(np.mean(_cell_fibers_densities))
            _offset_index += 1

    # pairs
    print('Pairs')
    _experiments = load.experiments_groups_as_tuples(PAIRS_EXPERIMENTS)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = filtering.by_real_cells(_experiments)
    _experiments = filtering.by_distance(_experiments, PAIRS_CELLS_DISTANCE)
    if PAIRS_BAND:
        _experiments = filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_x in OFFSETS_X:
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

    _rois = compute.rois(_arguments)
    _fibers_densities = compute.fibers_densities(_rois)

    _pairs_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Group:', _group, sep='\t')
        _offset_index = 0
        _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        # take offsets based on cells distance
        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_coordinates = [list(_properties['time_points'][0]['left_cell']['coordinates'].values())]
        _right_cell_coordinates = [list(_properties['time_points'][0]['right_cell']['coordinates'].values())]
        _cells_distance = compute.cells_distance_in_cell_size(
            _experiment, _series_id, _left_cell_coordinates, _right_cell_coordinates)
        _edges_distance = _cells_distance - 1
        _max_x_offset = _edges_distance - ROI_LENGTH

        for _offset_x in OFFSETS_X:
            if _offset_x > _max_x_offset:
                break

            _arguments = {
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
            }
            _roi = compute.roi_time_point(_arguments)
            _fibers_density = _fibers_densities[_roi]

            if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                continue

            _normalized_fibers_density = compute_lib.z_score(
                _x=_fibers_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )

            _pairs_fibers_densities[_offset_index].append(_normalized_fibers_density)
            _offset_index += 1

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[OFFSETS_X] * 2,
        _y_array=[_pairs_fibers_densities, _single_cell_fibers_densities],
        _names_array=['Pairs', 'Single Cell'],
        _modes_array=['lines+markers'] * 2,
        _dashes_array=['solid', 'dash'],
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score',
        _title='Fibers Densities vs. Distance from Cell'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_distance_' + str(PAIRS_CELLS_DISTANCE) + '_time_point_' + str(TIME_POINT)
    )


if __name__ == '__main__':
    main()
