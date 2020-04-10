import os
import sys

import numpy as np
import plotly.graph_objs as go
from scipy.stats import pearsonr

import libs.compute_lib
from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import save

EXPERIMENTS = ['SN16']
MINIMUM_TIME_POINTS = sys.maxsize
EARLY_TIME_POINTS = {
    True: {
        'start': {
            'SN16': 0,
            'SN18': 0,
            'SN41': 0,
            'SN44': 0
        },
        'end': {
            'SN16': 6,
            'SN18': 6,
            'SN41': 18,
            'SN44': 18
        }
    },
    False: {
        'start': {
            'SN16': 6,
            # TODO: update if needed
            'SN18': 0,
            'SN41': 0,
            'SN44': 0
        },
        'end': {
            'SN16': 12,
            # TODO: update if needed
            'SN18': 6,
            'SN41': 18,
            'SN44': 18
        }
    }
}
CELLS_DISTANCE_RANGE = (5, 18)
DIRECTION = 'inside'
REAL_CELLS = True
STATIC = False
BAND = True

OFFSET_X = 0.5
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 1

PLOT = True
CONDITIONAL_NORMALIZATION = False
X_LABELS_START = -3
X_LABELS_END = 14
Y_LABELS_START = -1.1
Y_LABELS_END = 3.4
X_BINS = 1
Y_BINS = 5
Z_MIN = 0
# according to conditional normalization
Z_MAX = {
    True: 0.2,
    False: 0.05
}


def main(_early_time_points=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
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
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _heatmap_fibers = []
    _heatmap_fibers_change = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _series_normalization = load.normalization_series_file_data(_experiment, 'Series ' + str(_series_id))
        _series_normalization = [_series_normalization['average'], _series_normalization['std']]
        for _cell_id in ['left_cell', 'right_cell']:
            _fibers_densities_by_time = [_fibers_densities[_tuple] for _tuple in
                                         _rois_dictionary[(_experiment, _series_id, _group, _cell_id)]]
            _cell_fibers_densities = \
                _fibers_densities_by_time[EARLY_TIME_POINTS[_early_time_points]['start'][_experiment]:
                                          EARLY_TIME_POINTS[_early_time_points]['end'][_experiment]]
            _properties = load.group_properties(_experiment, _series_id, _group)
            _cell_fibers_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids'][_cell_id], _cell_fibers_densities)
            _cell_fibers_densities = compute.longest_fibers_densities_ascending_sequence(_cell_fibers_densities)

            # fix if found nan
            if True in np.isnan(_cell_fibers_densities):
                _cell_fibers_densities = _cell_fibers_densities[:np.where(np.isnan(_cell_fibers_densities))[0][0]]

            # not enough data
            if len(_cell_fibers_densities) < DERIVATIVE + 1:
                continue

            _z_score_fibers_density = libs.compute_lib.z_score_fibers_densities_array(
                _cell_fibers_densities, _series_normalization
            )
            if _experiment in ['SN41', 'SN44']:
                for _start_index in [0, 1, 2]:
                    _heatmap_fibers += _z_score_fibers_density[_start_index::3][DERIVATIVE:]
                    _heatmap_fibers_change += compute_lib.derivative(
                        _z_score_fibers_density[_start_index::3], _n=DERIVATIVE
                    )
            else:
                _heatmap_fibers += _z_score_fibers_density[DERIVATIVE:]
                _heatmap_fibers_change += compute_lib.derivative(
                    _z_score_fibers_density, _n=DERIVATIVE
                )

    print(pearsonr(_heatmap_fibers, _heatmap_fibers_change))

    if PLOT:
        _y_shape = int(round((Y_LABELS_END - Y_LABELS_START) * Y_BINS))
        _x_shape = int(round((X_LABELS_END - X_LABELS_START) * X_BINS))
        _total_points = 0
        _z_array = np.zeros(shape=(_y_shape, _x_shape))
        for _y, _x in zip(_heatmap_fibers_change, _heatmap_fibers):
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
            for _fibers_index, _fibers_density_z_score in enumerate(_z_array):
                _sum = np.sum(_fibers_density_z_score)
                for _change_index, _change_z_score in enumerate(_fibers_density_z_score):
                    _z_array_plot[_fibers_index][_change_index] = (_change_z_score / _sum) if _sum != 0 else 0

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
                    'title': 'Fibers densities z-score',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Change in fibers<br>densities z-score',
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
                      str(BAND) + '_direction_' + DIRECTION + '_early_' + str(_early_time_points)
        )


if __name__ == '__main__':
    main()
