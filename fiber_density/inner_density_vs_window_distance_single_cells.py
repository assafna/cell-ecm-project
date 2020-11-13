import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib, paths_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import compute as experiments_compute
from libs.experiments import config as experiments_config
from libs.experiments import filtering as experiments_filtering
from libs.experiments import load as experiments_load
from libs.experiments import organize as experiments_organize
from libs.experiments.config import all_experiments, OUT_OF_BOUNDARIES
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering
from libs.simulations import load as simulations_load
from plotting import save

OFFSETS_X_END = 2.4
OFFSET_X_STEP = 0.2
OFFSETS_X = np.arange(start=0, stop=OFFSETS_X_END + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0

# experiments
EXPERIMENTS_TIME_FRAME = 12
OFFSET_Z = 0

# simulations
SIMULATIONS_TIME_POINT = {
    False: 50,
    True: 21
}


def compute_experiments_data():
    _experiments = all_experiments()
    _experiments = experiments_filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=True,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = experiments_load.experiments_groups_as_tuples(_experiments)
    _tuples = experiments_filtering.by_time_frames_amount(_tuples, EXPERIMENTS_TIME_FRAME)
    _tuples = experiments_filtering.by_main_cell(_tuples)

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right']):
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': experiments_config.QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': experiments_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': experiments_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': EXPERIMENTS_TIME_FRAME - 1
            })

    _windows_dictionary, _windows_to_compute = \
        experiments_compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'direction'])
    _fiber_densities = experiments_compute.fiber_densities(_windows_to_compute)

    _tuples = experiments_organize.by_single_cell_id(_tuples)

    _experiments_fiber_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in tqdm(_tuples, desc='Experiments loop'):
        _experiment, _series_id, _cell_id = _tuple
        _normalization = experiments_load.normalization_series_file_data(_experiment, _series_id)

        for _offset_x_index, _offset_x in enumerate(OFFSETS_X):
            _cell_fiber_densities = []
            for _cell_tuple in _tuples[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _offset_x, _direction)][0]
                    _fiber_density = _fiber_densities[_window_tuple]

                    if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                        continue

                    _normalized_fiber_density = compute_lib.z_score(
                        _x=_fiber_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    if not np.isnan(_normalized_fiber_density):
                        _cell_fiber_densities.append(_normalized_fiber_density)

            if len(_cell_fiber_densities) > 0:
                _experiments_fiber_densities[_offset_x_index].append(np.mean(_cell_fiber_densities))

    return _experiments_fiber_densities


def compute_simulations_fiber_densities(_simulations, _low_connectivity):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x, _direction in product(OFFSETS_X, ['left', 'right', 'up', 'down']):
            _arguments.append({
                'simulation': _simulation,
                'length_x': simulations_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else simulations_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': simulations_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else simulations_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else _offset_x,
                'offset_y': _offset_x if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_point': SIMULATIONS_TIME_POINT[_low_connectivity]
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.window_fiber_density_time_point, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            if _keys['direction'] in ['up', 'down']:
                _fiber_densities[(_keys['simulation'], _keys['offset_y'], _keys['direction'])] = _value
            else:
                _fiber_densities[(_keys['simulation'], _keys['offset_x'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def compute_simulations_data(_low_connectivity):
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(
        _simulations, _time_points=SIMULATIONS_TIME_POINT[_low_connectivity]
    )
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )

    _fiber_densities = compute_simulations_fiber_densities(_simulations, _low_connectivity)

    _simulations_fiber_densities = [[] for _i in range(len(OFFSETS_X))]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _normalization = simulations_load.normalization(_simulation)

        for _offset_x_index, _offset_x in enumerate(OFFSETS_X):
            _direction_fiber_densities = []
            for _direction in ['left', 'right', 'up', 'down']:
                _fiber_density = _fiber_densities[(_simulation, _offset_x, _direction)]

                _normalized_fiber_density = compute_lib.z_score(
                    _fiber_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fiber_densities.append(_normalized_fiber_density)

            _simulations_fiber_densities[_offset_x_index].append(np.mean(_direction_fiber_densities))

    return _simulations_fiber_densities


def main(_low_connectivity=False):
    print('Simulations')
    _simulations_fiber_densities = compute_simulations_data(_low_connectivity)

    print('Experiments')
    _experiments_fiber_densities = compute_experiments_data()

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _simulations_fiber_densities],
                name='Simulations',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _simulations_fiber_densities],
                    'thickness': 1,
                    'color': '#005b96'
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#005b96'
                },
                opacity=0.7
            ),
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _experiments_fiber_densities],
                name='Experiments',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _experiments_fiber_densities],
                    'thickness': 1,
                    'color': '#ea8500'
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': '#ea8500'
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-1.7, 13],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': 3.4,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': -0.2,
                    'y1': 13,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                }
            ]
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
