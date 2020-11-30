import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths, config
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, OUT_OF_BOUNDARIES, \
    all_experiments
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVES = [0, 1, 2]
DERIVATIVES_TEXT = ['D', 'D\'', 'D\'\'']


def compute_single_cell_mean(_experiment, _series_id, _cell_tuples, _windows_dictionary, _fiber_densities):
    _cell_fiber_densities = []
    for _time_frame in range(compute.density_time_frame(_experiment)):
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
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=True,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.density_time_frame(_experiments[0]))
    _tuples = filtering.by_main_cell(_tuples)

    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.density_time_frame(_experiment)
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
                'time_points': _time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'direction'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _tuples = organize.by_single_cell_id(_tuples)
    print('Total tuples:', len(_tuples))
    _experiments_ids = list(_tuples.keys())

    _y_arrays = [[] for _i in DERIVATIVES]
    for _index_1 in tqdm(range(len(_experiments_ids)), desc='Main loop'):
        _tuple_1 = _experiments_ids[_index_1]
        _experiment_1, _series_id_1, _cell_id_1 = _tuple_1
        _fiber_densities_1 = compute_single_cell_mean(
            _experiment=_experiment_1,
            _series_id=_series_id_1,
            _cell_tuples=_tuples[_tuple_1],
            _windows_dictionary=_windows_dictionary,
            _fiber_densities=_fiber_densities
        )
        for _index_2 in range(_index_1 + 1, len(_experiments_ids)):
            _tuple_2 = _experiments_ids[_index_2]
            _experiment_2, _series_id_2, _cell_id_2 = _tuple_2
            _fiber_densities_2 = compute_single_cell_mean(
                _experiment=_experiment_2,
                _series_id=_series_id_2,
                _cell_tuples=_tuples[_tuple_2],
                _windows_dictionary=_windows_dictionary,
                _fiber_densities=_fiber_densities
            )
            for _derivative_index, _derivative in enumerate(DERIVATIVES):
                _y_arrays[_derivative_index].append(compute_lib.correlation(
                    compute_lib.derivative(_fiber_densities_1, _n=_derivative),
                    compute_lib.derivative(_fiber_densities_2, _n=_derivative)
                ))

    print('Total points:', len(_y_arrays[0]))
    print('Wilcoxon around the zero')
    for _y_array, _derivative in zip(_y_arrays, DERIVATIVES):
        print('Derivative:', _derivative, wilcoxon(_y_array))

    # plot
    _colors_array = config.colors(3)
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
                'title': 'Correlation',
                'range': [-1, 1],
                'zeroline': False,
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot'
    )


if __name__ == '__main__':
    main()
