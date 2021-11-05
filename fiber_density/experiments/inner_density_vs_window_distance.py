import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

PAIR_DISTANCE_RANGES = [(4, 6), (6, 8), (8, 10)]
OFFSET_X_STEP = 0.2
OFFSET_Z = 0
OFFSET_Y = 0


def compute_cell_pairs():
    _x_array = []
    _y_array = []
    _names_array = []
    for _distances_range in PAIR_DISTANCE_RANGES:
        print('Pair distance range:', str(_distances_range))
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
        _tuples = filtering.by_pair_distance_range(_tuples, _distances_range)
        _tuples = filtering.by_band(_tuples)
        print('Total tuples:', len(_tuples))

        _max_offsets_x = []
        _arguments = []
        for _tuple in _tuples:
            _experiment, _series_id, _group = _tuple
            _time_frame = compute.minimum_time_frames_for_correlation(_experiment)
            _pair_distance = \
                compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame - 1)
            _offsets_x = np.arange(start=0,
                                   stop=_pair_distance / 2 - 0.5 - QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                                   step=OFFSET_X_STEP)
            if len(_offsets_x) > len(_max_offsets_x):
                _max_offsets_x = _offsets_x
            for _offset_x in _offsets_x:
                for _cell_id in ['left_cell', 'right_cell']:
                    _arguments.append({
                        'experiment': _experiment,
                        'series_id': _series_id,
                        'group': _group,
                        'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                        'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                        'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                        'offset_x': _offset_x,
                        'offset_y': OFFSET_Y,
                        'offset_z': OFFSET_Z,
                        'cell_id': _cell_id,
                        'direction': 'inside',
                        'time_point': _time_frame - 1
                    })

        _windows_dictionary, _windows_to_compute = \
            compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'offset_x', 'cell_id'])
        _fiber_densities = compute.fiber_densities(_windows_to_compute)

        _pair_distance_fiber_densities = [[] for _i in range(len(_max_offsets_x))]
        for _tuple in _tuples:
            _experiment, _series_id, _group = _tuple
            for _offset_x_index, _offset_x in enumerate(_max_offsets_x):
                for _cell_id in ['left_cell', 'right_cell']:
                    if (_experiment, _series_id, _group, _offset_x, _cell_id) in _windows_dictionary:
                        _normalization = load.normalization_series_file_data(_experiment, _series_id)
                        _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _offset_x, _cell_id)][0]
                        _fiber_density = _fiber_densities[_window_tuple]

                        if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                            continue

                        _normalized_fiber_density = compute_lib.z_score(
                            _x=_fiber_density[0],
                            _average=_normalization['average'],
                            _std=_normalization['std']
                        )

                        if not np.isnan(_normalized_fiber_density):
                            _pair_distance_fiber_densities[_offset_x_index].append(_normalized_fiber_density)

        _x_array.append(_max_offsets_x)
        _y_array.append(_pair_distance_fiber_densities)
        _names_array.append('Pair distance ' + str(_distances_range[0]) + '-' + str(_distances_range[1]))
    return _names_array, _x_array, _y_array


def main():
    _names_array, _x_array, _y_array = compute_cell_pairs()

    # plot
    _colors_array = config.colors(3)
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x,
                y=[np.mean(_array) for _array in _y],
                name=_name,
                error_y={
                    'type': 'data',
                    'array': [np.std(_array) for _array in _y],
                    'thickness': 1,
                    'color': _color
                },
                mode='markers',
                marker={
                    'size': 15,
                    'color': _color
                },
                opacity=0.7
            ) for _x, _y, _name, _color in zip(_x_array, _y_array, _names_array, _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'range': [-1.7, 13],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [0, 4, 8, 12]
            },
            'legend': {
                'xanchor': 'right',
                'yanchor': 'top',
                'bordercolor': 'black',
                'borderwidth': 2,
                'bgcolor': 'white'
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': 3.4,
                    'y1': -1.5,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -1.5,
                    'x1': -0.2,
                    'y1': 13,
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
        _filename='plot'
    )


if __name__ == '__main__':
    main()
