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
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering
from libs.simulations import load as simulations_load
from plotting import scatter, save, update

OFFSET_X = 0
OFFSET_Y = 0

# experiments
EXPERIMENTS_TIME_POINTS = 18
OFFSET_Z = 0
OUT_OF_BOUNDARIES = True

# simulations
SIMULATIONS_TIME_POINTS = 50


def compute_simulations_fibers_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _direction in ['left', 'right', 'up', 'down']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': simulations_config.ROI_HEIGHT
                if _direction in ['up', 'down'] else simulations_config.ROI_WIDTH,
                'length_y': simulations_config.ROI_WIDTH
                if _direction in ['up', 'down'] else simulations_config.ROI_HEIGHT,
                'offset_x': OFFSET_Y if _direction in ['up', 'down'] else OFFSET_X,
                'offset_y': OFFSET_X if _direction in ['up', 'down'] else OFFSET_Y,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': SIMULATIONS_TIME_POINTS
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.roi_fibers_density_by_time, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[(_keys['simulation'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main():
    # experiments
    print('Experiments')
    _experiments = experiments_load.experiments_groups_as_tuples(experiments_config.SINGLE_CELL)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINTS)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
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
                'cell_id': 'cell',
                'direction': _direction
            })

    _rois_dictionary, _rois_to_compute = \
        experiments_compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fibers_densities = experiments_compute.fibers_densities(_rois_to_compute)

    _experiments = experiments_organize.by_single_cell_id(_experiments)

    _experiments_fibers_densities = [[] for _i in range(EXPERIMENTS_TIME_POINTS)]
    for _tuple in tqdm(_experiments, desc='Experiments Loop'):
        _experiment, _series_id, _group = _tuple
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _time_point in range(EXPERIMENTS_TIME_POINTS):
            _cell_fibers_densities = []
            for _cell_tuple in _experiments[_tuple]:
                _, _, _group = _cell_tuple
                for _direction in ['left', 'right']:
                    _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _direction)][_time_point]
                    _fibers_density = _fibers_densities[_roi_tuple]

                    if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                        continue

                    _normalized_fibers_density = compute_lib.z_score(
                        _x=_fibers_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    if not np.isnan(_normalized_fibers_density):
                        _cell_fibers_densities.append(_normalized_fibers_density)

            if len(_cell_fibers_densities) > 0:
                _experiments_fibers_densities[_time_point].append(np.mean(_cell_fibers_densities))

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
        _is_dominant_passive=False
    )

    _fibers_densities = compute_simulations_fibers_densities(_simulations)

    _simulations_fibers_densities = [[] for _i in range(SIMULATIONS_TIME_POINTS)]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _normalization = simulations_load.normalization(_simulation)

        for _time_point in range(SIMULATIONS_TIME_POINTS):
            _direction_fibers_densities = []
            for _direction in ['left', 'right', 'up', 'down']:
                _fibers_density = _fibers_densities[(_simulation, _direction)][_time_point]

                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _direction_fibers_densities.append(_normalized_fibers_density)

            _simulations_fibers_densities[_time_point].append(np.mean(_direction_fibers_densities))

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=np.array(range(EXPERIMENTS_TIME_POINTS)) * 15,
                y=[np.mean(_array) for _array in _experiments_fibers_densities],
                name='Experiments',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _experiments_fibers_densities],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'dash'}
            ),
            go.Scatter(
                x=list(range(SIMULATIONS_TIME_POINTS)),
                xaxis='x2',
                y=[np.mean(_array) for _array in _simulations_fibers_densities],
                name='Simulations',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _simulations_fibers_densities],
                    'thickness': 1
                },
                mode='lines+markers',
                line={'dash': 'solid'}
            )
        ],
        layout={
            'xaxis': {
                'title': 'Time (minutes)',
                'titlefont': {
                    'color': 'rgb(31, 119, 180)'
                },
                'tickfont': {
                    'color': 'rgb(31, 119, 180)'
                },
                'side': 'top'
            },
            'xaxis2': {
                'title': 'Cell contraction (percentages)',
                'titlefont': {
                    'color': 'rgb(255, 127, 14)'
                },
                'tickfont': {
                    'color': 'rgb(255, 127, 14)'
                },
                # 'anchor': 'free',
                'overlaying': 'x',
                'side': 'bottom'
            },
            'yaxis': {
                'title': 'Fibers Density Z-score',
                # 'domain': [0.3, 1],
                'range': [-1.5, 9]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths_lib.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
