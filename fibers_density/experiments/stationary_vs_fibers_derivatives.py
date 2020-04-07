import os

import plotly.graph_objs as go
from statsmodels.tsa.stattools import kpss, adfuller
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH
from plotting import save

EXPERIMENTS = ['SN16']
CELLS_DISTANCE_RANGE = (6, 8)
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
BAND = True
OUT_OF_BOUNDARIES = False
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['I', 'I\'', 'I\'\'']


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments)
    if BAND:
        _experiments = filtering.by_band(_experiments)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _cell_id in ['left_cell', 'right_cell']:
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
                'cell_id': _cell_id,
                'direction': 'inside'
            })

    _rois_dictionary, _rois_to_compute = \
        compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _kpss_y_arrays = [[] for _i in DERIVATIVES]
    _adf_y_arrays = [[] for _i in DERIVATIVES]
    for _tuple in tqdm(_experiments, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple
        _properties = load.group_properties(_experiment, _series_id, _group)
        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fibers_densities = _experiments_fibers_densities[(_experiment, _series_id, _group, _cell_id)]
            _cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids'][_cell_id], _cell_fibers_densities)

            if not OUT_OF_BOUNDARIES:
                _cell_fibers_densities = compute.longest_fibers_densities_ascending_sequence(_cell_fibers_densities)
            else:
                _cell_fibers_densities = [_fibers_density[0] for _fibers_density in _cell_fibers_densities]

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
                        'color': _color,
                        # 'line': {
                        #     'width': 1,
                        #     'color': 'white'
                        # }
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
