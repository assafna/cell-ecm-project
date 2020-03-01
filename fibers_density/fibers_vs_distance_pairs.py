import os
from multiprocessing.pool import Pool

import numpy as np
from tqdm import tqdm

from libs import compute_lib, paths_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load as experiments_load
from libs.experiments import compute as experiments_compute
from libs.experiments import config as experiments_config
from libs.experiments import filtering as experiments_filtering

from libs.simulations import load as simulations_load
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering

from libs.simulations.config import CELL_DIAMETER
from plotting import scatter, save

CELLS_DISTANCE = 7
OFFSET_X_STEP = 0.1
OFFSET_Y = 0

# experiments
EXPERIMENTS = ['SN16']
EXPERIMENTS_TIME_POINT = 18
BAND = True
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False
EXPERIMENTS_OFFSET_X_END = {
    5: 2.8,
    7: 5.3
}
EXPERIMENTS_OFFSETS_X = \
    np.arange(start=0, stop=EXPERIMENTS_OFFSET_X_END[CELLS_DISTANCE] + OFFSET_X_STEP, step=OFFSET_X_STEP)

# simulations
SIMULATIONS_TIME_POINT = 50
SIMULATIONS_OFFSET_X_END = CELLS_DISTANCE - (SIMULATIONS_TIME_POINT / 100) - simulations_config.ROI_WIDTH
SIMULATIONS_OFFSETS_X = \
    np.arange(start=0, stop=SIMULATIONS_OFFSET_X_END + OFFSET_X_STEP, step=OFFSET_X_STEP)


def compute_experiments_fibers_densities(_experiments):
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_x in EXPERIMENTS_OFFSETS_X:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': experiments_config.ROI_LENGTH,
                'length_y': experiments_config.ROI_HEIGHT,
                'length_z': experiments_config.ROI_WIDTH,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'left_cell',
                'direction': 'inside',
                'time_point': EXPERIMENTS_TIME_POINT - 1,
                'save': False
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(experiments_compute.roi_fibers_density_time_point, _arguments), total=len(_arguments)):
            _fibers_densities[(_keys['experiment'], _keys['series_id'], _keys['group'], _keys['offset_x'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def compute_simulations_fibers_densities(_simulations):
    _fibers_densities = list([None] * len(SIMULATIONS_OFFSETS_X))
    for _simulation in _simulations:
        print(_simulation)
        _offset_index = 0
        _normalization = simulations_load.normalization(_simulation)
        for _offset_x in SIMULATIONS_OFFSETS_X:
            _fibers_density = simulations_compute.roi_fibers_density_time_point(
                _simulation=_simulation,
                _length_x=simulations_config.ROI_WIDTH,
                _length_y=simulations_config.ROI_HEIGHT,
                _offset_x=_offset_x,
                _offset_y=OFFSET_Y,
                _cell_id='left_cell',
                _direction='inside',
                _time_point=SIMULATIONS_TIME_POINT
            )
            _normalized_fibers_density = compute_lib.z_score(
                _fibers_density,
                _normalization['average'],
                _normalization['std']
            )

            if _fibers_densities[_offset_index] is None:
                _fibers_densities[_offset_index] = [_normalized_fibers_density]
            else:
                _fibers_densities[_offset_index].append(_normalized_fibers_density)
            _offset_index += 1

    return _fibers_densities


def main():
    # experiments
    _experiments = experiments_load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINT)
    _experiments = experiments_filtering.by_real_cells(_experiments)
    _experiments = experiments_filtering.by_distance(
        _experiments, CELLS_DISTANCE, _time_point=EXPERIMENTS_TIME_POINT - 1)
    if BAND:
        _experiments = experiments_filtering.by_band(_experiments)

    _fibers_densities = compute_experiments_fibers_densities(_experiments)

    _experiments_fibers_densities = [[] for _i in range(len(EXPERIMENTS_OFFSETS_X))]
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        print('Experiment:', _experiment, 'Series ID:', _series_id, 'Group:', _group, sep='\t')
        _offset_index = 0
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _offset_x in EXPERIMENTS_OFFSETS_X:
            _fibers_density = _fibers_densities[(_experiment, _series_id, _group, _offset_x)]
            if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                continue

            _normalized_fibers_density = compute_lib.z_score(
                _x=_fibers_density[0],
                _average=_normalization['average'],
                _std=_normalization['std']
            )

            _experiments_fibers_densities[_offset_index].append(_normalized_fibers_density)
            _offset_index += 1

    # simulations
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINT)
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = simulations_filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)

    _simulations_fibers_densities = compute_simulations_fibers_densities(_simulations)

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[np.arange(start=0, stop=1000, step=OFFSET_X_STEP)] * 2,
        _y_array=[_experiments_fibers_densities, _simulations_fibers_densities],
        _names_array=['Experiments', 'Simulations'],
        _modes_array=['lines+markers'] * 2,
        _dashes_array=['dash', 'solid'],
        _x_axis_title='Distance from Left Cell (cell size)',
        _y_axis_title='Fibers Density Z-score'
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
