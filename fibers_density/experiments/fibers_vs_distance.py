import os
from itertools import product

import numpy as np
from multiprocess.pool import Pool

from libs import compute_lib
from libs.config_lib import USE_MULTIPROCESSING, CPUS_TO_USE
from libs.experiments import load, compute, config, filtering, organize, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

PAIRS_EXPERIMENTS = ['SN16']
PAIRS_CELLS_DISTANCE = 7
PAIRS_BAND = True
TIME_POINT = 15
OFFSET_X_END = PAIRS_CELLS_DISTANCE - 1
OFFSET_X_STEP = 0.1
OFFSETS_X = np.arange(start=0, stop=OFFSET_X_END + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False


def compute_offsets(_experiment, _series_id, _group, _cell_id, _directions):
    for _offset_x, _direction in product(OFFSETS_X, _directions):
        compute.roi_fibers_density_time_point(_experiment, _series_id, _group, ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT,
                                              _offset_x, OFFSET_Y, OFFSET_Z, _cell_id, _direction, TIME_POINT - 1)


def main():
    # single cell
    print('Single Cell')
    _experiments = load.experiments_groups_as_tuples(config.SINGLE_CELL)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINT)
    _experiments = organize.by_single_cell_id(_experiments)

    # prepare data in mp
    if USE_MULTIPROCESSING:
        _arguments = []
        for _tuple in _experiments:
            _experiment, _series_id, _cell_id = _tuple
            for _cell_tuple in _experiments[_tuple]:
                _, _, _group = _cell_tuple
                _arguments.append((_experiment, _series_id, _group, 'cell', ['left', 'right']))

        _p = Pool(CPUS_TO_USE)
        _p.starmap(compute_offsets, _arguments)
        _p.close()

    _single_cell_fibers_densities = list([None] * len(OFFSETS_X))
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
                    _fibers_density = compute.roi_fibers_density_time_point(_experiment, _series_id, _group, ROI_LENGTH,
                                                                            ROI_WIDTH, ROI_HEIGHT, _offset_x, OFFSET_Y,
                                                                            OFFSET_Z, 'cell', _direction, TIME_POINT - 1
                                                                            )
                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    _cell_fibers_densities.append(_normalized_fibers_density)

            if len(_cell_fibers_densities) > 0:
                if _single_cell_fibers_densities[_offset_index] is None:
                    _single_cell_fibers_densities[_offset_index] = [np.mean(_cell_fibers_densities)]
                else:
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

    # prepare data in mp
    if USE_MULTIPROCESSING:
        _arguments = []
        for _tuple in _experiments:
            _experiment, _series_id, _group = _tuple
            _arguments.append((_experiment, _series_id, _group, 'left_cell', ['inside']))

        _p = Pool(CPUS_TO_USE)
        _p.starmap(compute_offsets, _arguments)
        _p.close()

    _pairs_fibers_densities = list([None] * len(OFFSETS_X))
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Group:', _group, sep='\t')
        _offset_index = 0
        _normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
        for _offset_x in OFFSETS_X:
            _fibers_density = compute.roi_fibers_density_time_point(_experiment, _series_id, _group, ROI_LENGTH,
                                                                    ROI_WIDTH, ROI_HEIGHT, _offset_x, OFFSET_Y,
                                                                    OFFSET_Z, 'left_cell', 'inside', TIME_POINT - 1
                                                                    )
            if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                continue

            _normalized_fibers_density = compute_lib.z_score(
                _x=_fibers_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )

            if _pairs_fibers_densities[_offset_index] is None:
                _pairs_fibers_densities[_offset_index] = [_normalized_fibers_density]
            else:
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
