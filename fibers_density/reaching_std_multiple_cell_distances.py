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
from libs.simulations import compute as simulations_compute
from libs.simulations import config as simulations_config
from libs.simulations import filtering as simulations_filtering
from libs.simulations import load as simulations_load
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0

# experiments
CELLS_DISTANCE_RANGES = [(4, 6), (6, 8), (8, 10)]
EXPERIMENTS = ['SN16']
EXPERIMENTS_TIME_POINTS = 18
BAND = True
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False
EXPERIMENTS_Z_SCORE_GOAL = 5

# simulations
CELLS_DISTANCES = [5, 7, 9]
SIMULATIONS_TIME_POINTS = 50
SIMULATIONS_STEP = int(round(SIMULATIONS_TIME_POINTS / EXPERIMENTS_TIME_POINTS))
SIMULATIONS_Z_SCORE_GOAL = 2


def compute_experiments_data(_cells_distance_range):
    _experiments = experiments_load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINTS)
    _experiments = experiments_filtering.by_real_cells(_experiments)
    _experiments = experiments_filtering.by_distance_range(_experiments, _cells_distance_range)
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
                'time_points': EXPERIMENTS_TIME_POINTS - 1
            })

    _rois_dictionary, _rois_to_compute = \
        experiments_compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = experiments_compute.fibers_densities(_rois_to_compute)

    _experiments_time_points = []
    for _tuple in tqdm(_experiments, desc='Experiments Loop'):
        _experiment, _series_id, _group = _tuple
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _time_point in range(EXPERIMENTS_TIME_POINTS - 1):
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
                    if _normalized_fibers_density > EXPERIMENTS_Z_SCORE_GOAL:
                        _experiments_time_points.append(max(0, _time_point - 1))
                        break
            else:
                continue
            break

    return _experiments_time_points


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


def compute_simulations_data(_cells_distance):
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
    _simulations = simulations_filtering.by_distance(_simulations, _distance=_cells_distance)

    _fibers_densities = compute_simulations_fibers_densities(_simulations)

    _simulations_fibers_densities = []
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

                if not np.isnan(_normalized_fibers_density):
                    if _normalized_fibers_density > SIMULATIONS_Z_SCORE_GOAL:
                        _simulations_fibers_densities.append(max(0, _time_point - 1))
                        break
            else:
                continue
            break

    return _simulations_fibers_densities


def main():
    print('Simulations')
    _simulations_y_array = []
    for _cells_distance in CELLS_DISTANCES:
        print('Cells Distance:', _cells_distance)
        _simulations_fibers_densities = compute_simulations_data(_cells_distance)
        _simulations_y_array.append(_simulations_fibers_densities)

    print('Experiments')
    _experiments_y_array = []
    for _cells_distance_range in CELLS_DISTANCE_RANGES:
        print('Cells Distance:', _cells_distance_range)
        _experiments_fibers_densities = compute_experiments_data(_cells_distance_range)
        _experiments_y_array.append(_experiments_fibers_densities)

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=np.array(CELLS_DISTANCES) - 0.1,
                y=[np.mean(_array) for _array in _simulations_y_array],
                name='Simulations',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _simulations_y_array],
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
                x=np.array(CELLS_DISTANCES) + 0.1,
                y=[np.mean(_array) * 15 for _array in _experiments_y_array],
                yaxis='y2',
                name='Experiments',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) * 15 for _array in _experiments_y_array],
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
                'title': 'Cells distance (cell diameter)',
                'tickmode': 'array',
                'tickvals': CELLS_DISTANCES
            },
            'yaxis': {
                'title': 'Cell contraction<br>to z-score ' + str(SIMULATIONS_Z_SCORE_GOAL) + ' (percentages)',
                'titlefont': {
                    'color': '#005b96'
                },
                'tickfont': {
                    'color': '#005b96'
                },
                'range': [14, 50],
                # 'anchor': 'free',
                'showgrid': False,
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [15, 25, 35, 45]
            },
            'yaxis2': {
                'title': 'Time to z-score ' + str(EXPERIMENTS_Z_SCORE_GOAL) + ' (minutes)',
                'titlefont': {
                    'color': '#ea8500'
                },
                'tickfont': {
                    'color': '#ea8500'
                },
                'range': [-20, 120],
                'overlaying': 'y',
                'side': 'right',
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 20, 60, 100]
            },
            'legend': {
                'xanchor': 'left',
                'x': 0.1,
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2
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
