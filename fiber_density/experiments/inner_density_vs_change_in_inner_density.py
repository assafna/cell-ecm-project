import os
import sys

import numpy as np
import plotly.graph_objs as go

import libs.compute_lib
from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER
from plotting import save

EXPERIMENTS = ['SN16']
MINIMUM_TIME_FRAMES = sys.maxsize
# when using offset 0.5
# TIME_FRAMES = {
#     'early': [0, 6],
#     'late': [17, 23]
# }
# when using offset 0:
TIME_FRAMES = {
    'early': [0, 6],
    'late': [12, 18]
}
PAIR_DISTANCE_RANGE = (5, 18)
DIRECTION = 'inside'
REAL_CELLS = True
STATIC = False
BAND = True

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 1

PLOT = True
CONDITIONAL_NORMALIZATION = False
X_LABELS_START = -3
X_LABELS_END = 14
Y_LABELS_START = -1.5
Y_LABELS_END = 3.4
X_BINS = 1
Y_BINS = 5
Z_MIN = 0
# according to conditional normalization
Z_MAX = {
    True: 0.2,
    False: 0.05
}


def main(_early_time_frames=True):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _experiments = filtering.by_real_pairs(_experiments, _real_pairs=REAL_CELLS)
    _experiments = filtering.by_fake_static_pairs(_experiments, _fake_static_pairs=STATIC)
    _experiments = filtering.by_band(_experiments, _band=BAND)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
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
                'direction': DIRECTION
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _heatmap_fiber = []
    _heatmap_fiber_change = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _series_normalization = load.normalization_series_file_data(_experiment, _series_id)
        for _cell_id in ['left_cell', 'right_cell']:
            _fiber_densities_by_time = [_fiber_densities[_tuple] for _tuple in
                                        _windows_dictionary[(_experiment, _series_id, _group, _cell_id)]]
            _cell_fiber_densities = \
                _fiber_densities_by_time[TIME_FRAMES['early'][0] if _early_time_frames else TIME_FRAMES['late'][0]:
                                         TIME_FRAMES['early'][1] if _early_time_frames else TIME_FRAMES['late'][1]]
            _properties = load.group_properties(_experiment, _series_id, _group)
            _cell_fiber_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids'][_cell_id], _cell_fiber_densities)
            _cell_fiber_densities = compute.longest_fiber_densities_ascending_sequence(_cell_fiber_densities)

            # fix if found nan
            if True in np.isnan(_cell_fiber_densities):
                _cell_fiber_densities = _cell_fiber_densities[:np.where(np.isnan(_cell_fiber_densities))[0][0]]

            # not enough data
            if len(_cell_fiber_densities) < DERIVATIVE + 1:
                continue

            _z_score_fiber_density = libs.compute_lib.z_score_array(
                _array=_cell_fiber_densities,
                _average=_series_normalization['average'],
                _std=_series_normalization['std']
            )
            if _experiment in ['SN41', 'SN44']:
                for _start_index in [0, 1, 2]:
                    _heatmap_fiber += _z_score_fiber_density[_start_index::3][DERIVATIVE:]
                    _heatmap_fiber_change += compute_lib.derivative(
                        _z_score_fiber_density[_start_index::3], _n=DERIVATIVE
                    )
            else:
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
                    'title': 'Fiber densities z-score',
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
            _filename='plot_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                      str(BAND) + '_direction_' + DIRECTION + '_early_' + str(_early_time_frames)
        )


if __name__ == '__main__':
    main()
