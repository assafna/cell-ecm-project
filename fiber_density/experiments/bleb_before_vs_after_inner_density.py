import os

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, \
    AFTER_BLEB_INJECTION_FIRST_TIME_FRAME, all_experiments
from plotting import save

OFFSET_X = 0
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False


def main(_band=True, _offset_y=0):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=None,
        _is_bleb=True,
        _is_bleb_from_start=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, _distance_range=PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=REAL_CELLS)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=STATIC)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

    _arguments = []
    _longest_time_frame = 0
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)

        # save for later
        if _latest_time_frame > _longest_time_frame:
            _longest_time_frame = _latest_time_frame

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': _offset_y,
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

    _valid_tuples = []
    _valid_cells = []
    _densities = [[] for _ in range(_longest_time_frame)]
    for _tuple in tqdm(_tuples, desc='Experiments loop'):
        _experiment, _series_id, _group = _tuple
        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        _properties = load.group_properties(_experiment, _series_id, _group)

        for _cell_id in ['left_cell', 'right_cell']:
            _cell_fiber_densities = \
                _experiments_fiber_densities[(_experiment, _series_id, _group, _cell_id)]

            # not enough time frames
            if len(_cell_fiber_densities) < AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiment]:
                continue

            _cell_fiber_densities = compute.remove_blacklist(
                _experiment, _series_id, _properties['cells_ids'][_cell_id], _cell_fiber_densities)

            for _time_frame, _cell_fiber_density in enumerate(_cell_fiber_densities):

                # not out of border
                if _cell_fiber_density[1]:
                    continue

                # normalize
                _cell_fiber_density_normalized = compute_lib.z_score(
                    _x=_cell_fiber_density[0],
                    _average=_normalization['average'],
                    _std=_normalization['std']
                )

                # save
                _densities[_time_frame].append(_cell_fiber_density_normalized)

                if _tuple not in _valid_tuples:
                    _valid_tuples.append(_tuple)

                _cell_tuple = (_experiment, _series_id, _group, _cell_id)
                if _cell_tuple not in _valid_cells:
                    _valid_cells.append(_cell_tuple)

    print('Total pairs:', len(_valid_tuples))
    print('Total cells:', len(_valid_cells))

    # plot
    _temporal_resolution = compute.temporal_resolution_in_minutes(_experiments[0])
    _first_time_frame_after_injection = AFTER_BLEB_INJECTION_FIRST_TIME_FRAME[_experiments[0]]
    _fig = go.Figure(
        data=go.Scatter(
            x=np.array(range(_longest_time_frame)) * _temporal_resolution,
            y=[np.mean(_array) for _array in _densities],
            name='Fiber density (z-score)',
            error_y={
                'type': 'data',
                'array': [np.std(_array) for _array in _densities],
                'thickness': 1
            },
            mode='lines+markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            line={'dash': 'solid'},
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Time (minutes)',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Fiber density (z-score)',
                'zeroline': False
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': -_temporal_resolution,
                    'y0': -1,
                    'x1': -_temporal_resolution,
                    'y1': 8,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': -_temporal_resolution,
                    'y0': -1,
                    'x1': 350,
                    'y1': -1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': (_first_time_frame_after_injection - 0.5) * _temporal_resolution,
                    'y0': -1,
                    'x1': (_first_time_frame_after_injection - 0.5) * _temporal_resolution,
                    'y1': 8,
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
        _filename='plot_band_' + str(_band) + '_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
