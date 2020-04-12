import os

import numpy as np
import plotly.graph_objs as go
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import SINGLE_CELL, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
TIME_POINTS = 18
OUT_OF_BOUNDARIES = False
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_single_cell_mean(_experiment, _series_id, _cell_tuples, _rois_dictionary, _fibers_densities):
    _cell_fibers_densities = []
    for _time_point in range(TIME_POINTS):
        _time_point_fibers_densities = []
        for _cell_tuple in _cell_tuples:
            _, _, _group = _cell_tuple
            for _direction in ['left', 'right']:
                _roi_tuple = _rois_dictionary[(_experiment, _series_id, _group, _direction)][_time_point]
                _fibers_density = _fibers_densities[_roi_tuple]

                if not OUT_OF_BOUNDARIES and _fibers_density[1]:
                    continue

                _time_point_fibers_densities.append(_fibers_density)

        _cell_fibers_densities.append(np.mean(_time_point_fibers_densities))

    return _cell_fibers_densities


def main():
    _experiments = load.experiments_groups_as_tuples(SINGLE_CELL)
    _experiments = filtering.by_time_points_amount(_experiments, TIME_POINTS)
    _experiments = filtering.by_main_cell(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _direction in ['left', 'right']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': 'cell',
                'direction': _direction,
                'time_points': TIME_POINTS
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments = organize.by_single_cell_id(_experiments)
    print('Total experiments:', len(_experiments))

    _kpss_y_arrays = [[] for _i in DERIVATIVES]
    _adf_y_arrays = [[] for _i in DERIVATIVES]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _cell_id = _tuple
        _cell_fibers_densities = compute_single_cell_mean(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_tuples=_experiments[_tuple],
            _rois_dictionary=_rois_dictionary,
            _fibers_densities=_fibers_densities
        )
        for _derivative_index, _derivative in enumerate(DERIVATIVES):
            _cell_fibers_densities_derivative = compute_lib.derivative(_cell_fibers_densities, _n=_derivative)
            _, _kpss_p_value, _, _ = kpss(_cell_fibers_densities_derivative)
            _kpss_y_arrays[_derivative_index].append(_kpss_p_value)
            _, _adf_p_value, _, _, _, _ = adfuller(_cell_fibers_densities_derivative)
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
        _stationary_count = len([_value for _value in _adf_y_arrays[_derivative_index] if _value < 0.95])
        print('Derivative:', _derivative, 'Stationary:',
              str(_stationary_count / len(_adf_y_arrays[_derivative_index]) * 100) + '%')

    # plot
    _colors_array = ['#844b00', '#ea8500', '#edbc80']
    for _test_name, _y_title, _y_tickvals, _p_value_line, _y_arrays in \
            zip(
                ['kpss', 'adf'],
                ['KPSS test p-value', 'ADF test p-value'],
                [[0.05, 0.1], [0, 0.5, 0.95]],
                [0.05, 0.95],
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
                    'title': 'Density by derivatives',
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
