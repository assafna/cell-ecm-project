import os
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save, update

TIME_POINT = {
    False: 50,
    True: 35
}
CELLS_DISTANCES = [5, 7, 9]
OFFSET_X_STEP = 0.2
OFFSET_X_END = {
    5: 1.6,
    7: 2.6,
    9: 3.6
}
OFFSET_Y = 0


def compute_simulations_fibers_densities(_simulations, _offsets_x, _low_connectivity):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x in _offsets_x:
            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'simulation': _simulation,
                    'length_x': ROI_WIDTH,
                    'length_y': ROI_HEIGHT,
                    'offset_x': _offset_x,
                    'offset_y': OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': TIME_POINT[_low_connectivity]
                })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_time_point, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[
                (_keys['simulation'], _keys['offset_x'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fibers_densities


def main(_low_connectivity=False):
    _x_array = []
    _y_array = []
    _names_array = []
    for _distance in CELLS_DISTANCES:
        print('Cells Distance ' + str(_distance))
        _simulations = load.structured()
        _simulations = filtering.by_categories(
            _simulations,
            _is_single_cell=False,
            _is_heterogeneity=False,
            _is_low_connectivity=_low_connectivity,
            _is_causality=False,
            _is_dominant_passive=False
        )
        _simulations = filtering.by_distance(_simulations, _distance=_distance)
        _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT[_low_connectivity])

        _offsets_x = np.arange(start=0, stop=OFFSET_X_END[_distance] + OFFSET_X_STEP, step=OFFSET_X_STEP)
        _fibers_densities = compute_simulations_fibers_densities(_simulations, _offsets_x, _low_connectivity)

        _cells_distance_fibers_densities = [[] for _i in range(len(_offsets_x))]
        for _simulation in _simulations:
            _offset_index = 0
            _normalization = load.normalization(_simulation)
            for _offset_x in _offsets_x:
                for _cell_id in ['left_cell', 'right_cell']:
                    _fibers_density = _fibers_densities[(_simulation, _offset_x, _cell_id)]

                    _normalized_fibers_density = compute_lib.z_score(
                        _fibers_density,
                        _normalization['average'],
                        _normalization['std']
                    )
                    _cells_distance_fibers_densities[_offset_index].append(_normalized_fibers_density)

                _offset_index += 1

        _x_array.append(_offsets_x)
        _y_array.append(_cells_distance_fibers_densities)
        _names_array.append('Distance ' + str(_distance))

    # plot
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=[np.mean(_array) for _array in _y],
                name=_name,
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _y],
                    'thickness': 2
                },
                mode='markers',
                marker={
                    'size': 15
                },
                opacity=0.7
            ) for _x, _y, _name in zip(_x_array, _y_array, _names_array)
        ],
        layout={
            'xaxis': {
                'title': 'Distance from Cell (cell size)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fibers Density Z-score',
                'range': [-1.7, 16],
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
                    'y0': -1.5,
                    'x1': max([OFFSET_X_END[_distance] for _distance in CELLS_DISTANCES]) + OFFSET_X_STEP,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -OFFSET_X_STEP,
                    'y0': -1.5,
                    'x1': -OFFSET_X_STEP,
                    'y1': 16,
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
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
