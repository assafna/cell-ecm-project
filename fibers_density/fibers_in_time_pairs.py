import os
from multiprocessing.pool import Pool

import numpy as np
from tqdm import tqdm

from libs import compute_lib, paths_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import compute as experiments_compute
from libs.experiments import config as experiments_config
from libs.experiments import filtering as experiments_filtering
from libs.experiments import load as experiments_load
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering
from libs.simulations import load as simulations_load
from plotting import scatter, save, edit

CELLS_DISTANCE = 7
OFFSET_X = 0
OFFSET_Y = 0

# experiments
EXPERIMENTS = ['SN16']
EXPERIMENTS_TIME_POINTS = 18
BAND = True
OFFSET_Z = 0
OUT_OF_BOUNDARIES = True

# simulations
SIMULATIONS_TIME_POINTS = 50


def compute_simulations_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': simulations_config.ROI_WIDTH,
                'length_y': simulations_config.ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': SIMULATIONS_TIME_POINTS
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[(_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    # experiments
    print('Experiments')
    _experiments = experiments_load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINTS)
    _experiments = experiments_filtering.by_real_cells(_experiments)
    _experiments = experiments_filtering.by_distance(_experiments, CELLS_DISTANCE)
    if BAND:
        _experiments = experiments_filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': experiments_config.ROI_LENGTH,
                'length_y': experiments_config.ROI_HEIGHT,
                'length_z': experiments_config.ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': EXPERIMENTS_TIME_POINTS
            })

    _rois_dictionary, _rois_to_compute = \
        experiments_compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = experiments_compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = [[] for _i in range(EXPERIMENTS_TIME_POINTS)]
    for _tuple in tqdm(_experiments, desc='Experiments Loop'):
        _experiment, _series_id, _group = _tuple
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _time_point in range(EXPERIMENTS_TIME_POINTS):
            for _cell_id in ['left_cell', 'right_cell']:
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _cell_id)][_time_point]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _normalized_fibers_density = compute_lib.z_score(
                    _x=_fibers_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                if not np.isnan(_normalized_fibers_density):
                    _experiments_fibers_densities[_time_point].append(_normalized_fibers_density)

    # simulations
    print('Simulations')
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINTS)
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = simulations_filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)

    _fibers_densities = compute_simulations_fibers_densities(_simulations)

    _simulations_fibers_densities = [[] for _i in range(SIMULATIONS_TIME_POINTS)]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _normalization = simulations_load.normalization(_simulation)

        for _time_point in range(SIMULATIONS_TIME_POINTS):
            for _cell_id in ['left_cell', 'right_cell']:
                _fibers_density = _fibers_densities[(_simulation, _cell_id)][_time_point]

                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )

                _simulations_fibers_densities[_time_point].append(_normalized_fibers_density)

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[
            np.arange(
                start=0,
                stop=SIMULATIONS_TIME_POINTS + 1,
                step=(SIMULATIONS_TIME_POINTS - 1) / (EXPERIMENTS_TIME_POINTS - 1)
            ),
            list(range(SIMULATIONS_TIME_POINTS))
        ],
        _y_array=[_experiments_fibers_densities, _simulations_fibers_densities],
        _names_array=['Experiments', 'Simulations'],
        _modes_array=['lines+markers'] * 2,
        _dashes_array=['dash', 'solid'],
        _x_axis_title='Time',
        _y_axis_title='Fibers Density Z-score'
    )

    _fig = edit.update_y_axis(
        _fig=_fig,
        _range=[-1.5, 9]
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
