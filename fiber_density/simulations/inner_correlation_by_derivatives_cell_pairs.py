import os
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths, config
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCE = 5
SIMULATIONS_TIME_POINTS = 50
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_fiber_densities(_simulations, _directions):
    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            for _direction in _directions:
                _arguments.append({
                    'simulation': _simulation,
                    'length_x': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                    'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                    'offset_x': OFFSET_X,
                    'offset_y': OFFSET_Y,
                    'cell_id': _cell_id,
                    'direction': _direction,
                    'time_points': SIMULATIONS_TIME_POINTS
                })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['cell_id'], _keys['direction'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main(_directions=None):
    if _directions is None:
        _directions = ['inside', 'outside']

    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=SIMULATIONS_TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_pair_distance(_simulations, _distance=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _fiber_densities = compute_fiber_densities(_simulations, _directions)

    for _direction in _directions:
        _y_arrays = [[] for _i in DERIVATIVES]
        for _simulation in tqdm(_simulations, desc='Simulations loop'):
            _left_cell_fiber_densities = _fiber_densities[(_simulation, 'left_cell', _direction)]
            _right_cell_fiber_densities = _fiber_densities[(_simulation, 'right_cell', _direction)]
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_left_cell_fiber_densities, _n=_derivative),
                    compute_lib.derivative(_right_cell_fiber_densities, _n=_derivative)
                ))

        print('Direction:', _direction)
        print('Total pairs:', len(_y_arrays[0]))
        print('Wilcoxon around the zero')
        for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
            print('Derivative:', _derivative, wilcoxon(_y_array))

        # plot
        _y_title = 'Inner correlation' if _direction == 'inside' else 'Outer correlation'
        _colors_array = config.colors(3)
        _fig = go.Figure(
            data=[
                go.Box(
                    y=_y,
                    name=_derivative,
                    boxpoints='all',
                    jitter=1,
                    pointpos=0,
                    line={
                        'width': 1
                    },
                    fillcolor='white',
                    marker={
                        'size': 10,
                        'color': _color
                    },
                    opacity=0.7,
                    showlegend=False
                ) for _y, _derivative, _color in zip(_y_arrays, DERIVATIVES_TEXT, _colors_array)
            ],
            layout={
                'xaxis': {
                    'title': 'Fiber density derivative',
                    'zeroline': False
                },
                'yaxis': {
                    'title': _y_title,
                    'range': [-1, 1],
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': [-1, -0.5, 0, 0.5, 1]
                }
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_direction_' + _direction
        )


if __name__ == '__main__':
    main()
