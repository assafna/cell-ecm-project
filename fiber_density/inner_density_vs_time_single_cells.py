import os
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

OFFSET_X = 0
OFFSET_Y = 0

# experiments
EXPERIMENTS_TIME_FRAMES = 12
EXPERIMENTS_TEMPORAL_RESOLUTION = 21
OFFSET_Z = 0

# simulations
SIMULATIONS_TIME_POINTS = 50
SIMULATIONS_STEP = int(round(SIMULATIONS_TIME_POINTS / EXPERIMENTS_TIME_FRAMES))


def compute_simulations_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _direction in ['left', 'right', 'up', 'down']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': simulations_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else simulations_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': simulations_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER
                if _direction in ['up', 'down'] else simulations_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                'offset_y': OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': SIMULATIONS_TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    # simulations
    print('Simulations')
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINTS)
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=True,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_simulations_fiber_densities(_simulations)

    _simulations_fiber_densities = [[] for _i in range(SIMULATIONS_TIME_POINTS)]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _normalization = simulations_load.normalization(_simulation)

        for _time_frame in range(SIMULATIONS_TIME_POINTS):
            _direction_fiber_densities = []
            for _direction in ['left', 'right', 'up', 'down']:
                _fiber_density = _fiber_densities[(_simulation, _direction)][_time_frame]

                _normalized_fiber_density = compute_lib.z_score(
                    _fiber_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fiber_densities.append(_normalized_fiber_density)

            _simulations_fiber_densities[_time_frame].append(np.mean(_direction_fiber_densities))

    print('Total simulations cells:', len(_simulations_fiber_densities[0]))

    # experiments
    print('Experiments')
    _experiments = all_experiments()
    _experiments = experiments_filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=True,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = experiments_load.experiments_groups_as_tuples(_experiments)
    _tuples = experiments_filtering.by_time_frames_amount(_tuples, EXPERIMENTS_TIME_FRAMES)
    _tuples = experiments_filtering.by_main_cell(_tuples)

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': experiments_config.QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': experiments_config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': experiments_config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': EXPERIMENTS_TIME_FRAMES
            })

    _windows_dictionary, _windows_to_compute = \
        experiments_compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fiber_densities = experiments_compute.fiber_densities(_windows_to_compute)

    _tuples = experiments_organize.by_single_cell_id(_tuples)
    print('Total experiments:', len(_tuples))

    _experiments_fiber_densities = [[] for _i in range(EXPERIMENTS_TIME_FRAMES)]
    for _tuple in tqdm(_tuples, desc='Experiments loop'):
        _experiment, _series_id, _ = _tuple
        _normalization = experiments_load.normalization_series_file_data(_experiment, _series_id)

        for _time_frame in range(EXPERIMENTS_TIME_FRAMES):
            _cell_fiber_densities = []
            for _cell_tuple in _tuples[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _direction)][_time_frame]
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
                _experiments_fiber_densities[_time_frame].append(np.mean(_cell_fiber_densities))

    print('Total experiments cells:', len(_experiments_fiber_densities[0]))

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=list(range(SIMULATIONS_TIME_POINTS))[::SIMULATIONS_STEP],
                y=[np.mean(_array) for _array in _simulations_fiber_densities][::SIMULATIONS_STEP],
                name='Simulations',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _simulations_fiber_densities][::SIMULATIONS_STEP],
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
                x=np.array(range(EXPERIMENTS_TIME_FRAMES)) * EXPERIMENTS_TEMPORAL_RESOLUTION,
                xaxis='x2',
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
                'title': 'Cell contraction (%)',
                'titlefont': {
                    'color': '#005b96'
                },
                'tickfont': {
                    'color': '#005b96'
                },
                'zeroline': False
            },
            'xaxis2': {
                'title': 'Time (minutes)',
                'titlefont': {
                    'color': '#ea8500'
                },
                'tickfont': {
                    'color': '#ea8500'
                },
                # 'anchor': 'free',
                'overlaying': 'x',
                'side': 'bottom',
                'showgrid': False,
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                # 'domain': [0.3, 1],
                'range': [-1.7, 14],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'left',
                'x': 0.1,
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -2,
                    'y0': -1.5,
                    'x1': 53,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -2,
                    'y0': -1.5,
                    'x1': -2,
                    'y1': 14,
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
