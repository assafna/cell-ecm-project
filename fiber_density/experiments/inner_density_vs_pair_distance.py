import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, compute, filtering, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

OFFSETS_X = [1, 2.6]
OFFSET_Y = 0
OFFSET_Z = 0


def main():
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.density_time_frame(_experiments[0]))
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.density_time_frame(_experiment)
        _pair_distance = \
            compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame=_time_frame - 1)
        for _offset_x in OFFSETS_X:
            if _pair_distance / 2 - 0.5 - QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER >= _offset_x:
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
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'offset_x'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    for _offset_x in OFFSETS_X:
        _x_array = []
        _y_array = []
        for _tuple in _tuples:
            _experiment, _series_id, _group = _tuple
            for _cell_id in ['left_cell', 'right_cell']:
                if (_experiment, _series_id, _group, _cell_id, _offset_x) in _windows_dictionary:
                    _pair_distance = \
                        compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame=0)
                    _normalization = load.normalization_series_file_data(_experiment, _series_id)
                    _window_tuple = _windows_dictionary[(_experiment, _series_id, _group, _cell_id, _offset_x)][0]
                    _fiber_density = _fiber_densities[_window_tuple]

                    if not OUT_OF_BOUNDARIES and _fiber_density[1]:
                        continue

                    _normalized_fiber_density = compute_lib.z_score(
                        _x=_fiber_density[0],
                        _average=_normalization['average'],
                        _std=_normalization['std']
                    )

                    if not np.isnan(_normalized_fiber_density):
                        _x_array.append(_pair_distance)
                        _y_array.append(_normalized_fiber_density)

        print('Offset x (cell diameter):', _offset_x)
        print('Total pairs:', len(_x_array))
        print(compute_lib.correlation(_x_array, _y_array, _with_p_value=True))

        # plot
        _fig = go.Figure(
            data=go.Scatter(
                x=_x_array,
                y=_y_array,
                mode='markers',
                marker={
                    'size': 15,
                    'color': 'black'
                }
            ),
            layout={
                'xaxis': {
                    'title': 'Pair distance (cell diameter)',
                    'zeroline': False
                },
                'yaxis': {
                    'title': 'Fiber density (z-score)',
                    'zeroline': False,
                    'range': [-2.2, 13],
                    'tickmode': 'array',
                    'tickvals': [0, 4, 8, 12]
                },
                'shapes': [
                    {
                        'type': 'line',
                        'x0': 4.5,
                        'y0': -2,
                        'x1': 9.5,
                        'y1': -2,
                        'line': {
                            'color': 'black',
                            'width': 2
                        }
                    },
                    {
                        'type': 'line',
                        'x0': 4.5,
                        'y0': -2,
                        'x1': 4.5,
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
            _filename='plot_offset_x_' + str(_offset_x)
        )


if __name__ == '__main__':
    main()
