import os
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

TIME_POINT = {
    False: 50,
    True: 35
}
PAIR_DISTANCE = [5, 7, 9]
OFFSET_X_STEP = 0.2
OFFSET_Y = 0


def compute_fiber_densities(_simulations, _offsets_x, _low_connectivity):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x in _offsets_x:
            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'simulation': _simulation,
                    'length_x': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                    'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                    'offset_x': _offset_x,
                    'offset_y': OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_point': TIME_POINT[_low_connectivity]
                })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_time_point, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[
                (_keys['simulation'], _keys['offset_x'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main(_low_connectivity=False):
    _x_array = []
    _y_array = []
    _names_array = []
    _max_offsets_x = []
    for _distance in PAIR_DISTANCE:
        print('Pair distance ' + str(_distance))
        _simulations = load.structured()
        _simulations = filtering.by_categories(
            _simulations,
            _is_single_cell=False,
            _is_heterogeneity=False,
            _is_low_connectivity=_low_connectivity,
            _is_causality=False,
            _is_dominant_passive=False
        )
        _simulations = filtering.by_pair_distance(_simulations, _distance=_distance)
        _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT[_low_connectivity])

        _offsets_x = np.arange(start=0, stop=_distance / 2 - 0.25 - QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                               step=OFFSET_X_STEP)
        if len(_offsets_x) > len(_max_offsets_x):
            _max_offsets_x = _offsets_x
        _fiber_densities = compute_fiber_densities(_simulations, _offsets_x, _low_connectivity)

        _pair_distance_fiber_densities = [[] for _i in range(len(_offsets_x))]
        for _simulation in _simulations:
            _offset_index = 0
            _normalization = load.normalization(_simulation)
            for _offset_x in _offsets_x:
                for _cell_id in ['left_cell', 'right_cell']:
                    _fiber_density = _fiber_densities[(_simulation, _offset_x, _cell_id)]

                    _normalized_fiber_density = compute_lib.z_score(
                        _fiber_density,
                        _normalization['average'],
                        _normalization['std']
                    )
                    _pair_distance_fiber_densities[_offset_index].append(_normalized_fiber_density)

                _offset_index += 1

        _x_array.append(_offsets_x)
        _y_array.append(_pair_distance_fiber_densities)
        _names_array.append('Pair distance ' + str(_distance))

    # plot
    _colors_array = ['#011f4b', '#005b96', '#74c2e8']
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=[np.mean(_array) for _array in _y],
                name=_name,
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _y],
                    'thickness': 1,
                    'color': _color
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _color
                },
                opacity=0.7
            ) for _x, _y, _name, _color in zip(_x_array, _y_array, _names_array, _colors_array)
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
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_low_con_' + str(_low_connectivity)
    )


if __name__ == '__main__':
    main()
