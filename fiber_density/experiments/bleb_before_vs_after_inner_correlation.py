import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, \
    AFTER_BLEB_INJECTION_FIRST_TIME_FRAME, DERIVATIVE, all_experiments
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0.5
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False


def main(_band=True):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=None,
        _is_bleb=True,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_experiments, _distance_range=PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=REAL_CELLS)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=STATIC)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in _tuples:
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

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _n_pairs = 0
    _before_correlations = []
    _after_correlations = []
    for _tuple in tqdm(_tuples, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple

        _left_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fiber_densities = \
            _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        _before_left_cell_fiber_densities = \
            _left_cell_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiment]]
        _before_right_cell_fiber_densities = \
            _right_cell_fiber_densities[:AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiment]]

        _after_left_cell_fiber_densities = \
            _left_cell_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiment]:]
        _after_right_cell_fiber_densities = \
            _right_cell_fiber_densities[AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiment]:]

        _before_left_cell_fiber_densities_filtered, _before_right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _before_left_cell_fiber_densities, _before_right_cell_fiber_densities)

        _after_left_cell_fiber_densities_filtered, _after_right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _after_left_cell_fiber_densities, _after_right_cell_fiber_densities)

        # ignore small arrays
        _minimum_time_frame_for_correlation = compute.minimum_time_frames_for_correlation(_experiment)
        if len(_before_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation or \
                len(_after_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation:
            continue

        _n_pairs += 1

        _before_correlations.append(compute_lib.correlation(
            compute_lib.derivative(_before_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_before_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        ))
        _after_correlations.append(compute_lib.correlation(
            compute_lib.derivative(_after_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_after_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        ))

    print('Total pairs:', _n_pairs)
    _before_minus_after = np.array(_before_correlations) - np.array(_after_correlations)
    print('Wilcoxon of before minus after around the zero:')
    print(wilcoxon(_before_minus_after))
    print('Higher before amount:', (_before_minus_after > 0).sum() / len(_before_minus_after))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_before_correlations,
            y=_after_correlations,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Correlation before bleb',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Correlation after bleb',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': -1,
                    'y1': 1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': 1,
                    'y1': -1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -1,
                    'y0': -1,
                    'x1': 1,
                    'y1': 1,
                    'line': {
                        'color': 'red',
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
