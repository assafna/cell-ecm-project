import os

import numpy as np
import plotly.graph_objs as go
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import SINGLE_CELL, QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
TIME_FRAMES = 18
OUT_OF_BOUNDARIES = False
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_single_cell_mean(_experiment, _series_id, _cell_tuples, _windows_dictionary, _fiber_densities):
    _cell_fiber_densities = []
    for _time_frame in range(TIME_FRAMES):
        _time_frame_fiber_densities = []
        for _cell_tuple in _cell_tuples:
            _, _, _group = _cell_tuple
            for _direction in ['left', 'right']:
                _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _direction)][_time_frame]
                _fiber_density = _fiber_densities[_window_tuple]

                if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                    continue

                _time_frame_fiber_densities.append(_fiber_density)

        _cell_fiber_densities.append(np.mean(_time_frame_fiber_densities))

    return _cell_fiber_densities


def main():
    _experiments = load.experiments_groups_as_tuples(SINGLE_CELL)
    _experiments = filtering.by_time_frames_amount(_experiments, TIME_FRAMES)
    _experiments = filtering.by_main_cell(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_FRAMES
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments = organize.by_single_cell_id(_experiments)
    print('Total experiments:', len(_experiments))

    _kpss_y_arrays = [[] for _i in DERIVATIVES]
    _adf_y_arrays = [[] for _i in DERIVATIVES]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _cell_id = _tuple
        _cell_fiber_densities = compute_single_cell_mean(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_tuples=_experiments[_tuple],
            _windows_dictionary=_windows_dictionary,
            _fiber_densities=_fiber_densities
        )
        for _derivative_index, _derivative in enumerate(DERIVATIVES):
            _cell_fiber_densities_derivative = compute_lib.derivative(_cell_fiber_densities, _n=_derivative)
            _, _kpss_p_value, _, _ = kpss(_cell_fiber_densities_derivative, nlags='legacy')
            _kpss_y_arrays[_derivative_index].append(_kpss_p_value)
            _, _adf_p_value, _, _, _, _ = adfuller(_cell_fiber_densities_derivative)
            _adf_y_arrays[_derivative_index].append(_adf_p_value)

    print('Total cells:', len(_kpss_y_arrays[0]))

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
    _colors_array = ['#844b00', '#ea8500', '#edbc80']
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
