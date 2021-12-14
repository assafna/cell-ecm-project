import os

import numpy as np
import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Z = 0

PAIR_DISTANCE_RANGE = [4, 10]


def main(_real_cells=True, _static=False, _dead_dead=False, _live_dead=False, _dead=False, _live=False, _bead=False,
         _metastasis=False, _bleb=False, _bleb_amount_um=None, _band=True, _offset_y=0.5,
         _high_temporal_resolution=False):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=_bleb,
        _is_dead_dead=_dead_dead,
        _is_live_dead=_live_dead,
        _is_bead=_bead,
        _is_metastasis=_metastasis
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=_real_cells)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=_static)
    if _dead_dead is not False or _live_dead is not False:
        _tuples = filtering.by_dead_live(_tuples, _dead=_dead, _live=_live)
    _tuples = filtering.by_band(_tuples, _band=_band)
    if _bleb:
        _tuples = filtering.by_bleb_amount_um(_tuples, _amount_um=_bleb_amount_um)
    _tuples_matched = organize.by_matched_real_real_and_real_fake(_tuples)
    print('Total matched pairs:', len(_tuples_matched))

    _arguments = []
    for _matched_tuple in _tuples_matched:
        for _tuple in _matched_tuple:
            _experiment, _series_id, _group = _tuple
            _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, OFFSET_X)
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
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _real_real_correlations_array = []
    _real_fake_correlations_array = []
    _valid_real_real_tuples = []
    for _tuple_matched in tqdm(_tuples_matched, desc='Main loop'):
        _real_real_tuple, _real_fake_tuple = _tuple_matched

        _real_real_experiment, _real_real_series, _real_real_group = _real_real_tuple
        _, _, _real_fake_group = _real_fake_tuple

        _real_real_left_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_real_real_experiment, _real_real_series, _real_real_group, 'left_cell')
            ]
        _real_real_right_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_real_real_experiment, _real_real_series, _real_real_group, 'right_cell')
            ]
        _real_fake_left_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_real_real_experiment, _real_real_series, _real_fake_group, 'left_cell')
            ]
        _real_fake_right_cell_fiber_densities = \
            _experiments_fiber_densities[
                (_real_real_experiment, _real_real_series, _real_fake_group, 'right_cell')
            ]

        _properties = \
            load.group_properties(_real_real_experiment, _real_real_series, _real_real_group)

        _real_real_left_cell_fiber_densities = compute.remove_blacklist(
            _real_real_experiment,
            _real_real_series,
            _properties['cells_ids']['left_cell'],
            _real_real_left_cell_fiber_densities
        )
        _real_real_right_cell_fiber_densities = compute.remove_blacklist(
            _real_real_experiment,
            _real_real_series,
            _properties['cells_ids']['right_cell'],
            _real_real_right_cell_fiber_densities
        )
        _real_fake_left_cell_fiber_densities = compute.remove_blacklist(
            _real_real_experiment,
            _real_real_series,
            _properties['cells_ids']['left_cell'],
            _real_fake_left_cell_fiber_densities
        )
        _real_fake_right_cell_fiber_densities = compute.remove_blacklist(
            _real_real_experiment,
            _real_real_series,
            _properties['cells_ids']['right_cell'],
            _real_fake_right_cell_fiber_densities
        )

        _real_real_left_cell_fiber_densities_filtered, _real_real_right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _real_real_left_cell_fiber_densities, _real_real_right_cell_fiber_densities
            )

        # ignore small arrays
        _minimum_time_frame_for_correlation = compute.minimum_time_frames_for_correlation(_real_real_experiment)
        if len(_real_real_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation:
            continue

        _real_real_correlation = compute_lib.correlation(
            compute_lib.derivative(_real_real_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_real_real_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        )

        _real_fake_left_cell_fiber_densities_filtered, _real_fake_right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _real_fake_left_cell_fiber_densities, _real_fake_right_cell_fiber_densities
            )

        # ignore small arrays
        if len(_real_fake_left_cell_fiber_densities_filtered) < _minimum_time_frame_for_correlation:
            continue

        _real_fake_correlation = compute_lib.correlation(
            compute_lib.derivative(_real_fake_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_real_fake_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        )

        _real_real_correlations_array.append(_real_real_correlation)
        _real_fake_correlations_array.append(_real_fake_correlation)

        if _real_real_tuple not in _valid_real_real_tuples:
            _valid_real_real_tuples.append(_real_real_tuple)

    print('Total real-real pairs:', len(_valid_real_real_tuples))
    _real_minus_fake = np.array(_real_real_correlations_array) - np.array(_real_fake_correlations_array)
    print('Wilcoxon of real-real minus real-fake around the zero:')
    print(wilcoxon(_real_minus_fake))
    print('Higher real-real amount:', (_real_minus_fake > 0).sum() / len(_real_minus_fake))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_real_real_correlations_array,
            y=_real_fake_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Real-real correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Real-fake correlation',
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
        _filename='plot_real_' + str(_real_cells) + '_static_' + str(_static) + '_dead_dead_' + str(_dead_dead) +
                  '_live_dead_' + str(_live_dead) + '_dead_' + str(_dead) + '_live_' + str(_live) + '_bead_' +
                  str(_bead) + '_metastasis_' + str(_metastasis) + '_bleb_' + str(_bleb) + str(_bleb_amount_um) +
                  '_band_' + str(_band) + '_high_time_' + str(_high_temporal_resolution) + '_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
