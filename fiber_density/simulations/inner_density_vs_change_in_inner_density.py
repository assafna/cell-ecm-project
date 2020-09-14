import os
import sys
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import filtering, load, compute, paths, organize
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

MINIMUM_TIME_POINTS = sys.maxsize
PAIR_DISTANCE = range(5, 18)
DIRECTION = 'inside'
HETEROGENEITY = False

OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1

PLOT = True
CONDITIONAL_NORMALIZATION = False
X_LABELS_START = -4
X_LABELS_END = 10.5
Y_LABELS_START = -0.4
Y_LABELS_END = 0.65
X_BINS = 2
Y_BINS = 20
Z_MIN = 0
# according to conditional normalization
Z_MAX = {
    True: 0.2,
    False: 0.05
}


def compute_fiber_densities(_simulations):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'cell_id': _cell_id,
                'direction': DIRECTION,
                'time_points': MINIMUM_TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[
                (_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=HETEROGENEITY,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_pair_distances(_simulations, PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))
    _simulations_by_distances = organize.by_pair_distance(_simulations)
    for _distance in _simulations_by_distances:
        print('Distance ', _distance, ', total simulations:', len(_simulations_by_distances[_distance]))

    _fiber_densities = compute_fiber_densities(_simulations)

    _heatmap_fiber = []
    _heatmap_fiber_change = []
    for _simulation in _simulations:
        _simulation_normalization = load.normalization(_simulation)
        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fiber_densities = _fiber_densities[(_simulation, _cell_id)]

            # not enough data
            if len(_cell_fiber_densities) < DERIVATIVE + 1:
                continue

            _z_score_fiber_density = compute_lib.z_score_array(
                _array=_cell_fiber_densities,
                _average=_simulation_normalization['average'],
                _std=_simulation_normalization['std']
            )
            _heatmap_fiber += _z_score_fiber_density[DERIVATIVE:]
            _heatmap_fiber_change += compute_lib.derivative(
                _z_score_fiber_density, _n=DERIVATIVE
            )

    print(compute_lib.correlation(_heatmap_fiber, _heatmap_fiber_change, _with_p_value=True))

    if PLOT:
        _y_shape = int(round((Y_LABELS_END - Y_LABELS_START) * Y_BINS))
        _x_shape = int(round((X_LABELS_END - X_LABELS_START) * X_BINS))
        _total_points = 0
        _z_array = np.zeros(shape=(_y_shape, _x_shape))
        for _y, _x in zip(_heatmap_fiber_change, _heatmap_fiber):
            _y_rounded, _x_rounded = int(round(_y * Y_BINS)), int(round(_x * X_BINS))
            _y_index, _x_index = int(_y_rounded - Y_LABELS_START * Y_BINS), int(_x_rounded - X_LABELS_START * X_BINS)
            if 0 <= _y_index < _z_array.shape[0] and 0 <= _x_index < _z_array.shape[1]:
                _z_array[_y_index][_x_index] += 1
                _total_points += 1
        _z_array = _z_array / _total_points

        if not CONDITIONAL_NORMALIZATION:
            _z_array[_z_array == 0] = None
        else:
            _z_array_plot = np.zeros(shape=np.array(_z_array).shape)
            for _fiber_index, _fiber_density_z_score in enumerate(_z_array):
                _sum = np.sum(_fiber_density_z_score)
                for _change_index, _change_z_score in enumerate(_fiber_density_z_score):
                    _z_array_plot[_fiber_index][_change_index] = (_change_z_score / _sum) if _sum != 0 else 0

            _z_array_plot[_z_array_plot == 0] = None

        _fig = go.Figure(
            data=go.Heatmap(
                x=np.arange(start=X_LABELS_START, stop=X_LABELS_END, step=1 / X_BINS),
                y=np.arange(start=Y_LABELS_START, stop=Y_LABELS_END, step=1 / Y_BINS),
                z=_z_array,
                colorscale='Viridis',
                colorbar={
                    'tickmode': 'array',
                    'tickvals': [0, 0.025, 0.05],
                    'ticktext': ['0', 'Fraction', '0.05'],
                    'tickangle': -90
                },
                zmin=Z_MIN,
                zmax=Z_MAX[CONDITIONAL_NORMALIZATION]
            ),
            layout={
                'xaxis': {
                    'title': 'fiber densities z-score',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Change in fiber<br>density (z-score)',
                    'zeroline': False
                },
                'shapes': [
                    {
                        'type': 'line',
                        'x0': X_LABELS_START,
                        'y0': Y_LABELS_START,
                        'x1': X_LABELS_END,
                        'y1': Y_LABELS_START,
                        'line': {
                            'color': 'black',
                            'width': 2
                        }
                    },
                    {
                        'type': 'line',
                        'x0': X_LABELS_START,
                        'y0': Y_LABELS_START,
                        'x1': X_LABELS_START,
                        'y1': Y_LABELS_END,
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
            _filename='plot_direction_' + DIRECTION
        )


if __name__ == '__main__':
    main()
