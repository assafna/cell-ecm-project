import os

import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

PAIR_DISTANCE_RANGE = [4, 10]

OFFSET_X = 0
OFFSET_Y_DENSITY = 0
OFFSET_Y_CORRELATION = 0.5
OFFSET_Z = 0


def main(_real_cells=True, _static=False, _band=True, _high_temporal_resolution=False):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=_high_temporal_resolution,
        _is_bleb=False,
        _is_bleb_from_start=False,
        _is_dead_live=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=_real_cells)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=_static)
    _tuples = filtering.by_band(_tuples, _band=_band)
    _tuples = filtering.by_time_frames_amount(_tuples, compute.minimum_time_frames_for_correlation(_experiments[0]))
    print('Total tuples:', len(_tuples))

    print('Density')
    _arguments = []
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _time_frame = compute.minimum_time_frames_for_correlation(_experiment)
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                'length_y': QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'length_z': QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y_DENSITY,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_point': _time_frame - 1
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _densities_fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=False)

    _densities_experiments_fiber_densities = {
        _key: [_densities_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    print('Correlations')
    _arguments = []
    for _tuple in _tuples:
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
                'offset_y': OFFSET_Y_CORRELATION,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _correlations_fiber_densities = compute.fiber_densities(_windows_to_compute, _subtract_border=True)

    _correlations_experiments_fiber_densities = {
        _key: [_correlations_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _densities = []
    _correlations = []
    for _tuple in tqdm(_tuples, desc='Main loop'):
        _experiment, _series_id, _group = _tuple

        # density
        _left_cell_fiber_density = \
            _densities_experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')][0]
        _right_cell_fiber_density = \
            _densities_experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')][0]

        # not relevant
        if _left_cell_fiber_density[1] or _right_cell_fiber_density[1]:
            continue

        _normalization = load.normalization_series_file_data(_experiment, _series_id)
        _left_cell_fiber_density_normalized = compute_lib.z_score(
            _x=_left_cell_fiber_density[0],
            _average=_normalization['average'],
            _std=_normalization['std']
        )
        _right_cell_fiber_density_normalized = compute_lib.z_score(
            _x=_right_cell_fiber_density[0],
            _average=_normalization['average'],
            _std=_normalization['std']
        )

        _density = (_left_cell_fiber_density_normalized + _right_cell_fiber_density_normalized) / 2

        # correlation
        _left_cell_fiber_densities = \
            _correlations_experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell')]
        _right_cell_fiber_densities = \
            _correlations_experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell')]

        _properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['left_cell'], _left_cell_fiber_densities)
        _right_cell_fiber_densities = compute.remove_blacklist(
            _experiment, _series_id, _properties['cells_ids']['right_cell'], _right_cell_fiber_densities)

        _left_cell_fiber_densities_filtered, _right_cell_fiber_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _left_cell_fiber_densities, _right_cell_fiber_densities
            )

        # ignore small arrays
        if len(_left_cell_fiber_densities_filtered) < compute.minimum_time_frames_for_correlation(_experiment):
            continue

        _correlation = compute_lib.correlation(
            compute_lib.derivative(_left_cell_fiber_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_fiber_densities_filtered, _n=DERIVATIVE)
        )

        _densities.append(_density)
        _correlations.append(_correlation)

    print('Total tuples:', len(_densities))
    print('Pearson correlation of densities and correlations:')
    print(compute_lib.correlation(_densities, _correlations, _with_p_value=True))

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_densities,
            y=_correlations,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Fiber density (z-score)',
                'zeroline': False,
                'range': [-1.1, 15.2],
                # 'tickmode': 'array',
                # 'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
            'shapes': [
                {
                    'type': 'line',
                    'x0': 0,
                    'y0': -1,
                    'x1': 0,
                    'y1': 1,
                    'line': {
                        'color': 'black',
                        'width': 2
                    }
                },
                {
                    'type': 'line',
                    'x0': 0,
                    'y0': -1,
                    'x1': 15,
                    'y1': -1,
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
        _filename='plot_real_' + str(_real_cells) + '_static_' + str(_static) + '_band_' + str(_band) +
                  '_high_time_' + str(_high_temporal_resolution) + '_y_density_' + str(OFFSET_Y_DENSITY) +
                  '_y_correlation_' + str(OFFSET_Y_CORRELATION)
    )


if __name__ == '__main__':
    main()
