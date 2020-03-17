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

CELLS_DISTANCE = 7
OFFSET_X_STEP = 0.2
OFFSET_X_END = {
    5: 1.6,
    7: 2.6,
    9: 3.6
}
OFFSETS_X = np.arange(start=0, stop=OFFSET_X_END[CELLS_DISTANCE] + OFFSET_X_STEP, step=OFFSET_X_STEP)
OFFSET_Y = 0

# experiments
EXPERIMENTS = ['SN16']
EXPERIMENTS_TIME_POINT = 18
BAND = True
OFFSET_Z = 0
OUT_OF_BOUNDARIES = False

# simulations
SIMULATIONS_TIME_POINT = {
    False: 50,
    True: 35
}


def compute_experiments_data():
    _experiments = experiments_load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = experiments_filtering.by_time_points_amount(_experiments, EXPERIMENTS_TIME_POINT)
    _experiments = experiments_filtering.by_real_cells(_experiments)
    _experiments = experiments_filtering.by_distance(
        _experiments, CELLS_DISTANCE, _time_point=EXPERIMENTS_TIME_POINT - 1)
    if BAND:
        _experiments = experiments_filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _offset_x in OFFSETS_X:
            for _cell_id in ['left_cell', 'right_cell']:
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
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': EXPERIMENTS_TIME_POINT - 1
                })

    _rois_dictionary, _rois_to_compute = \
        experiments_compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'cell_id'])
    _fibers_densities = experiments_compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _tuple in tqdm(_experiments, desc='Experiments Loop'):
        _experiment, _series_id, _group = _tuple
        _offset_index = 0
        _normalization = experiments_load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))

        for _offset_x in OFFSETS_X:
            for _cell_id in ['left_cell', 'right_cell']:
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _offset_x, _cell_id)][0]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _normalized_fibers_density = compute_lib.z_score(
                    _x=_fibers_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )
                _experiments_fibers_densities[_offset_index].append(_normalized_fibers_density)

            _offset_index += 1

    return _experiments_fibers_densities


def compute_simulations_fibers_densities(_simulations, _low_connectivity):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x in OFFSETS_X:
            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'simulation': _simulation,
                    'length_x': simulations_config.ROI_WIDTH,
                    'length_y': simulations_config.ROI_HEIGHT,
                    'offset_x': _offset_x,
                    'offset_y': OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': SIMULATIONS_TIME_POINT[_low_connectivity]
                })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(simulations_compute.roi_fibers_density_time_point, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[
                (_keys['simulation'], _keys['offset_x'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def compute_simulations_data(_low_connectivity):
    _simulations = simulations_load.structured()
    _simulations = simulations_filtering.by_time_points_amount(_simulations, SIMULATIONS_TIME_POINT[_low_connectivity])
    _simulations = simulations_filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = simulations_filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)

    _fibers_densities = compute_simulations_fibers_densities(_simulations, _low_connectivity)

    _simulations_fibers_densities = [[] for _i in range(len(OFFSETS_X))]
    for _simulation in tqdm(_simulations, desc='Simulations Loop'):
        _offset_index = 0
        _normalization = simulations_load.normalization(_simulation)

        for _offset_x in OFFSETS_X:
            for _cell_id in ['left_cell', 'right_cell']:
                _fibers_density = _fibers_densities[(_simulation, _offset_x, _cell_id)]

                _normalized_fibers_density = compute_lib.z_score(
                    _fibers_density,
                    _normalization['average'],
                    _normalization['std']
                )
                _simulations_fibers_densities[_offset_index].append(_normalized_fibers_density)

            _offset_index += 1

    return _simulations_fibers_densities


def main(_low_connectivity=False):
    print('Simulations')
    _simulations_fibers_densities = compute_simulations_data(_low_connectivity)

    print('Experiments')
    _experiments_fibers_densities = compute_experiments_data()

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _simulations_fibers_densities],
                name='Simulations',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _simulations_fibers_densities],
                    'thickness': 1
                },
                mode='markers',
                marker={
                    'size': 15
                },
                opacity=0.7
            ),
            go.Scatter(
                x=OFFSETS_X,
                y=[np.mean(_array) for _array in _experiments_fibers_densities],
                name='Experiments',
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _experiments_fibers_densities],
                    'thickness': 1
                },
                mode='markers',
                marker={
                    'size': 15
                },
                opacity=0.7
            )
        ],
        layout={
            'xaxis': {
                'title': 'Distance from Cell (cell size)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fibers Density Z-score',
                'range': [-1.2, 10],
                'zeroline': False
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
                    'x0': -OFFSET_X_STEP,
                    'y0': -1,
                    'x1': OFFSETS_X[-1] + OFFSET_X_STEP,
                    'y1': -1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -OFFSET_X_STEP,
                    'y0': -1,
                    'x1': -OFFSET_X_STEP,
                    'y1': 10,
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
        _filename='plot_cells_distance_' + str(CELLS_DISTANCE) + '_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
