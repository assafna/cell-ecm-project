import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import save

# by low connectivity
TIME_POINT = {
    False: 50,
    True: 35
}
OFFSET_X = 0
OFFSET_Y_START = -10
OFFSET_Y_END = 10
OFFSET_Y_STEP = 0.5
DERIVATIVE = 2
# by low connectivity
STD = {
    False: 0.5,
    True: 0.25
}
CELLS_DISTANCE = 5

# globals
_simulations = []
_fibers_densities = {}
_array = None


def compute_data(_arguments):
    _offset_y_index, _offset_y = _arguments

    _fibers_densities_array = []
    for _simulation in _simulations:
        _simulation_normalization = load.normalization(_simulation)
        for _cell_id in ['left_cell', 'right_cell']:
            _fibers_density = _fibers_densities[(_simulation, _offset_y, _cell_id)]
            _normalized_fibers_density = compute_lib.z_score(
                _x=_fibers_density,
                _average=_simulation_normalization['average'],
                _std=_simulation_normalization['std']
            )
            _fibers_densities_array.append(_normalized_fibers_density)

    if len(_fibers_densities_array) > 0:
        return _offset_y_index, np.mean(_fibers_densities_array)
    else:
        return _offset_y_index, None


def compute_array(_low_connectivity):
    global _simulations, _fibers_densities, _array

    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, TIME_POINT[_low_connectivity])
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=_low_connectivity,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD[_low_connectivity])
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    print('Total simulations:', len(_simulations))

    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _arguments = []
    for _simulation in _simulations:
        for _offset_y, _cell_id in product(_offsets_y, ['left_cell', 'right_cell']):
            _arguments.append({
                'simulation': _simulation,
                'length_x': ROI_WIDTH,
                'length_y': ROI_HEIGHT,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_point': TIME_POINT[_low_connectivity]
            })

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.roi_fibers_density_time_point, _arguments),
                total=len(_arguments), desc='Computing Rois & Fibers Densities'):
            _fibers_densities[(_keys['simulation'], _keys['offset_y'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    _arguments = []
    for _offset_y_index, _offset_y in enumerate(_offsets_y):
        _arguments.append((_offset_y_index, _offset_y))

    _array = np.zeros(shape=len(_offsets_y))
    with Pool(CPUS_TO_USE) as _p:
        for _answer in tqdm(_p.imap_unordered(compute_data, _arguments), total=len(_arguments),
                            desc='Computing array'):
            _offset_y_index, _mean = _answer
            _array[_offset_y_index] = _mean
        _p.close()
        _p.join()

    return _array


def main(_low_connectivity=False):
    global _simulations, _fibers_densities, _array

    compute_array(_low_connectivity)

    # plot
    _offsets_y = np.arange(start=OFFSET_Y_START, stop=OFFSET_Y_END + OFFSET_Y_STEP, step=OFFSET_Y_STEP)
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_offsets_y,
                y=_array,
                mode='markers+lines',
                marker={
                    'size': 10,
                    'color': 'black'
                },
                line={
                    'width': 1
                },
                showlegend=False
            )
        ],
        layout={
            'xaxis': {
                'title': 'Offset in Y axis (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score) mean',
                'zeroline': False,
                'range': [-1.1, 7],
                'tickmode': 'array',
                'tickvals': [0, 3, 6]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': OFFSET_Y_START,
                    'y0': -1,
                    'x1': OFFSET_Y_END,
                    'y1': -1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': OFFSET_Y_START,
                    'y0': -1,
                    'x1': OFFSET_Y_START,
                    'y1': 7,
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
