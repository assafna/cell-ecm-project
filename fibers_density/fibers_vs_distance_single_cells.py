import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from tqdm import tqdm

from libs import compute_lib, paths_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import compute as experiments_compute
from libs.experiments import config as experiments_config
from libs.experiments import filtering as experiments_filtering
from libs.experiments import load as experiments_load
from libs.experiments import organize as experiments_organize
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering
from libs.simulations import load as simulations_load
from plotting import scatter, save, edit

OFFSET_X_STEP = 0.1
OFFSET_X_END = 7
OFFSETS_X = np.arange(start=0, stop=OFFSET_X_END + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0

# experiments
EXPERIMENTS_TIME_POINT = 18
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False

# simulations
SIMULATIONS_TIME_POINT = 50


def compute_experiments_fibers_densities(_experiments):
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right']):
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
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': EXPERIMENTS_TIME_POINT - 1
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(experiments_compute.roi_fibers_density_time_point, _arguments), total=len(_arguments)):
            _fibers_densities[(_keys['experiment'], _keys['series_id'], _keys['group'], _keys['offset_x'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def compute_simulations_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        _offset_index = 0
        _normalization = simulations_load.normalization(_simulation)
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right', 'up', 'down']):
            _arguments.append({
                'simulation': _simulation,
                'length_x': simulations_config.ROI_HEIGHT
                if _direction in ['up', 'down'] else simulations_config.ROI_WIDTH,
                'length_y': simulations_config.ROI_WIDTH
                if _direction in ['up', 'down'] else simulations_config.ROI_HEIGHT,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else _offset_x,
                'offset_y': _offset_x if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': SIMULATIONS_TIME_POINT
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.roi_fibers_density_time_point, _arguments),
                total=len(_arguments)):
            if _keys['direction'] in ['up', 'down']:
                _fibers_densities[(_keys['simulation'], _keys['offset_y'], _keys['direction'])] = _value
            else:
                _fibers_densities[(_keys['simulation'], _keys['offset_x'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    # experiments
    print('Experiments')
    _experiments = experiments_load.experiments_groups_as_tuples(experiments_config.SINGLE_CELL)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINT)

    _fibers_densities = compute_experiments_fibers_densities(_experiments)

    _experiments = experiments_organize.by_single_cell_id(_experiments)

    _experiments_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _offset_index = 0
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _offset_x in OFFSETS_X:
            _cell_fibers_densities = []
            for _cell_tuple in _experiments[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _fibers_density = _fibers_densities[(_experiment, _series_id, _group, _offset_x, _direction)]
                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    _cell_fibers_densities.append(_normalized_fibers_density)

            if len(_cell_fibers_densities) > 0:
                _experiments_fibers_densities[_offset_index].append(np.mean(_cell_fibers_densities))
            _offset_index += 1

    # simulations
    print('Simulations')
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINT)
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )

    _fibers_densities = compute_simulations_fibers_densities(_simulations)

    _simulations_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _simulation in _simulations:
        _offset_index = 0
        _normalization = simulations_load.normalization(_simulation)

        for _offset_x in OFFSETS_X:
            _direction_fibers_densities = []
            for _direction in ['left', 'right', 'up', 'down']:
                _fibers_density = _fibers_densities[(_simulation, _offset_x, _direction)]

                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fibers_densities.append(_normalized_fibers_density)

            _simulations_fibers_densities[_offset_index].append(np.mean(_direction_fibers_densities))
            _offset_index += 1

    # plot
    _fig = scatter.create_error_bars_plot(
        _x_array=[OFFSETS_X] * 2,
        _y_array=[_experiments_fibers_densities, _simulations_fibers_densities],
        _names_array=['Experiments', 'Simulations'],
        _modes_array=['lines+markers'] * 2,
        _dashes_array=['dash', 'solid'],
        _x_axis_title='Distance from Cell (cell size)',
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
