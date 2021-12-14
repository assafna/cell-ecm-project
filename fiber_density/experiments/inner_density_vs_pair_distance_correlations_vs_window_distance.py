import os

import numpy as np
import plotly.graph_objs as go

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, all_experiments, \
    OUT_OF_BOUNDARIES
from plotting import save

OFFSET_X_STEP = 0.2
OFFSET_Z = 0
OFFSET_Y = 0


def main():
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.density_time_frame(_experiments[0]))
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)

    _max_offsets_x = []
    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.minimum_time_frames_for_correlation(_experiment)
        _pair_distance = \
            compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame=_time_frame - 1)
        _offsets_x = \
            np.arange(start=0, stop=_pair_distance / 2 - 0.5 - QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
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

    _x_array = []
    _y_array = []
    _n_array = []
    _p_value_array = []
    for _offset_x in _max_offsets_x:
        _pair_distances = []
        _z_scores = []
        for _tuple in _tuples:
            _experiment, _series_id, _group = _tuple
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
                        _pair_distances.append(
                            compute.pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group,
                                                                          _time_frame=0)
                        )
                        _z_scores.append(_normalized_fiber_density)

        if len(_pair_distances) > 2:
            _x_array.append(round(_offset_x, 1))
            _correlation = compute_lib.correlation(_pair_distances, _z_scores, _with_p_value=True)
            _y_array.append(round(_correlation[0], 2))
            _n_array.append(len(_pair_distances))
            _p_value_array.append(round(_correlation[1], 2))

    print('Pearson:')
    print(compute_lib.correlation(_x_array, _y_array, _with_p_value=True))

    # plot
    _significant_x_array = [_x for _x, _p_value in zip(_x_array, _p_value_array) if _p_value < 0.05]
    _fig = go.Figure(
        data=[
            go.Scatter(
                x=_x_array,
                y=_y_array,
                mode='markers',
                marker={
                    'size': 15,
                    'color': 'black'
                },
                showlegend=False
            ),
            go.Scatter(
                x=_significant_x_array,
                y=[-0.79] * len(_significant_x_array),
                mode='text',
                text='*',
                textfont={
                    'color': 'red'
                },
                showlegend=False
            )
        ],
        layout={
            'xaxis': {
                'title': 'Window distance (cell diameter)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Correlation:<br>pair distance vs. fiber density',
                'zeroline': False,
                'range': [-0.82, 0.3],
                'tickmode': 'array',
                'tickvals': [-0.75, -0.25, 0.25]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -0.8,
                    'x1': 3.2,
                    'y1': -0.8,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -0.2,
                    'y0': -0.8,
                    'x1': -0.2,
                    'y1': 0.3,
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

    # table
    print('Window distance (cell diameter)', 'Correlation: pair distance vs. fiber density', 'N', 'P-value',
          sep='\t')
    for _x, _y, _n, _p_value in zip(_x_array, _y_array, _n_array, _p_value_array):
        print(_x, _y, _n, _p_value, sep='\t')


if __name__ == '__main__':
    main()
