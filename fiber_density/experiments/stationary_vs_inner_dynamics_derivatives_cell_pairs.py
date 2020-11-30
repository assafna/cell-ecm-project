import os

import plotly.graph_objs as go
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]

DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def main():
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.density_time_frame(_experiments[0]))
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
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
    for _tuple in tqdm(_tuples, desc='Experiments loop'):
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
        _minimum_time_frames_for_correlation = compute.minimum_time_frames_for_correlation(_experiment)
        if len(_left_cell_fiber_densities) < _minimum_time_frames_for_correlation or \
                len(_right_cell_fiber_densities) < _minimum_time_frames_for_correlation:
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
