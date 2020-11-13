import os
import warnings
from multiprocessing.pool import Pool

import plotly.graph_objs as go
from statsmodels.tools.sm_exceptions import InterpolationWarning
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, compute, paths, config
from libs.simulations.config import QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
PAIR_DISTANCE = 7
SIMULATIONS_TIME_POINTS = 50
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


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
                'direction': 'inside',
                'time_points': SIMULATIONS_TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
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

    _fiber_densities = compute_fiber_densities(_simulations)

    _kpss_y_arrays = [[] for _i in DERIVATIVES]
    _adf_y_arrays = [[] for _i in DERIVATIVES]
    for _simulation in tqdm(_simulations, desc='Simulations loop'):
        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fiber_densities = _fiber_densities[(_simulation, _cell_id)]
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _cell_fiber_densities_derivative = compute_lib.derivative(_cell_fiber_densities, _n=_derivative)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', category=InterpolationWarning)
                    _, _kpss_p_value, _, _ = kpss(_cell_fiber_densities_derivative, nlags='legacy')
                    _kpss_y_arrays[_derivative_index].append(_kpss_p_value)
                    _, _adf_p_value, _, _, _, _ = adfuller(_cell_fiber_densities_derivative)
                    _adf_y_arrays[_derivative_index].append(_adf_p_value)

    # print results
    print('KPSS:')
    for _derivative_index, _derivative in enumerate(DERIVATIVES):
        _stationary_count = len([_value for _value in _kpss_y_arrays[_derivative_index] if _value > 0.05])
        print('Derivative:', _derivative, 'Stationary:',
              str(_stationary_count / len(_kpss_y_arrays[_derivative_index]) * 100) + '%')
    print('ADF:')
    for _derivative_index, _derivative in enumerate(DERIVATIVES):
        _stationary_count = len([_value for _value in _adf_y_arrays[_derivative_index] if _value < 0.05])
        print('Derivative:', _derivative, 'Stationary:',
              str(_stationary_count / len(_adf_y_arrays[_derivative_index]) * 100) + '%')

    # plot
    _colors_array = config.colors(3)
    for _test_name, _y_title, _y_tickvals, _p_value_line, _y_arrays in \
            zip(
                ['kpss', 'adf'],
                ['KPSS test p-value', 'ADF test p-value'],
                [[0.05, 0.1], [0.05, 1]],
                [0.05, 0.05],
                [_kpss_y_arrays, _adf_y_arrays]
            ):
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
                    'zeroline': False,
                    'tickmode': 'array',
                    'tickvals': _y_tickvals
                },
                'shapes': [
                    {
                        'type': 'line',
                        'x0': DERIVATIVES[0] - 0.75,
                        'y0': _p_value_line,
                        'x1': DERIVATIVES[-1] + 0.75,
                        'y1': _p_value_line,
                        'line': {
                            'color': 'red',
                            'width': 2,
                            'dash': 'dash'
                        }
                    }
                ]
            }
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_' + _test_name
        )


if __name__ == '__main__':
    main()
