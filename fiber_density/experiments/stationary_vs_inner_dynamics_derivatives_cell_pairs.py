import os

import plotly.graph_objs as go
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER
from plotting import save

EXPERIMENTS = ['SN16']
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
BAND = True
OUT_OF_BOUNDARIES = False
PAIR_DISTANCE_RANGE = [4, 10]
TIME_FRAMES = 18
MINIMUM_TIME_FRAMES = 15
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_time_frames_amount(_experiments, TIME_FRAMES)
    _experiments = filtering.by_real_pairs(_experiments)
    _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_frame = len(_properties['time_points'])
        for _time_frame in range(len(_properties['time_points'])):
            _pair_distance = \
                compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame)
            if _pair_distance - 1 - OFFSET_X * 2 < QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER * 2:
                _latest_time_frame = _time_frame - 1
                break

        for _cell_id in ['left_cell', 'right_cell']:
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
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _kpss_y_arrays = [[] for _i in DERIVATIVES]
    _adf_y_arrays = [[] for _i in DERIVATIVES]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple
        _properties = load.group_properties(_experiment, _series_id, _group)

        _left_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')]
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')]
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        if not OUT_OF_BOUNDARIES:
            _left_cell_fiber_densities = \
                compute.longest_fiber_densities_ascending_sequence(_left_cell_fiber_densities)
            _right_cell_fiber_densities = \
                compute.longest_fiber_densities_ascending_sequence(_right_cell_fiber_densities)
        else:
            _left_cell_fiber_densities = [_fiber_density[0] for _fiber_density in _left_cell_fiber_densities]
            _right_cell_fiber_densities = [_fiber_density[0] for _fiber_density in _right_cell_fiber_densities]

        # ignore small arrays
        if len(_left_cell_fiber_densities) < MINIMUM_TIME_FRAMES or len(_right_cell_fiber_densities) < \
                MINIMUM_TIME_FRAMES:
            continue

        for _derivative_index, _derivative in enumerate(DERIVATIVES):
            for _cell_fiber_densities in [_left_cell_fiber_densities, _right_cell_fiber_densities]:
                _cell_fiber_densities_derivative = compute_lib.derivative(_cell_fiber_densities, _n=_derivative)
                _, _kpss_p_value, _, _ = kpss(_cell_fiber_densities_derivative, nlags='legacy')
                _kpss_y_arrays[_derivative_index].append(_kpss_p_value)
                _, _adf_p_value, _, _, _, _ = adfuller(_cell_fiber_densities_derivative)
                _adf_y_arrays[_derivative_index].append(_adf_p_value)

    print('Total pairs:', len(_kpss_y_arrays[0]) / 2)

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
