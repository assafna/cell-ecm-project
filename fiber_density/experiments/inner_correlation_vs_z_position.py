import os

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths, organize
from libs.experiments.config import QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER, \
    QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER, QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER, all_experiments, \
    DERIVATIVE
from plotting import save

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]


def compute_fiber_densities(_band=True, _high_temporal_resolution=False):
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
    _tuples = filtering.by_pair_distance_range(_experiments, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples, _band=_band)
    print('Total tuples:', len(_tuples))

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

    _tuples_by_experiment = organize.by_experiment(_tuples)

    _correlations_array = []
    _z_positions_array = []
    for _experiment in _tuples_by_experiment:
        print('Experiment:', _experiment)
        _experiment_tuples = _tuples_by_experiment[_experiment]

        for _tuple in tqdm(_experiment_tuples, desc='Main loop'):
            _, _series, _group = _tuple

            _left_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series, _group, 'left_cell')]
            _right_cell_fiber_densities = _experiments_fiber_densities[(_experiment, _series, _group, 'right_cell')]

            _properties = load.group_properties(_experiment, _series, _group)
            _left_cell_fiber_densities = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['left_cell'],
                _left_cell_fiber_densities
            )
            _right_cell_fiber_densities = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['right_cell'],
                _right_cell_fiber_densities
            )

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

            _group_mean_z_position = compute.group_mean_z_position_from_substrate(_experiment, _series, _group)

            _correlations_array.append(_correlation)
            _z_positions_array.append(_group_mean_z_position)

    print('Total points:', len(_correlations_array))
    print('Wilcoxon of correlations around the zero:')
    print(wilcoxon(_correlations_array))
    print('Pearson correlation of correlations and z position distances:')
    print(compute_lib.correlation(_correlations_array, _z_positions_array, _with_p_value=True))

    return _correlations_array, _z_positions_array


def main(_band=True, _high_temporal_resolution=False):
    _correlations_array, _z_positions_array = compute_fiber_densities(_band, _high_temporal_resolution)

    _min_z_position, _max_z_position = min(_z_positions_array), max(_z_positions_array)
    _z_positions_array_normalized = [(_z_position_value - _min_z_position) / (_max_z_position - _min_z_position)
                                     for _z_position_value in _z_positions_array]

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_z_positions_array_normalized,
            y=_correlations_array,
            mode='markers',
            marker={
                'size': 5,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'Normalized mean Z distance',
                'zeroline': False,
                'range': [-0.1, 1.2],
                'tickmode': 'array',
                'tickvals': [0, 0.5, 1]
            },
            'yaxis': {
                'title': 'Inner window correlation',
                'zeroline': False,
                'range': [-1.1, 1.2],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            },
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_band_' + str(_band) + '_high_temporal_res_' + str(_high_temporal_resolution)
    )


if __name__ == '__main__':
    main()
