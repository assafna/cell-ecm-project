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
OFFSET_Z = 0
PAIR_DISTANCE_RANGE = [4, 10]


def compute_fiber_densities(_offset_y):
    _experiments = all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=False,
        _is_bleb=False,
        _is_bleb_from_start=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples)
    _tuples = filtering.by_band(_tuples)
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
                'offset_y': _offset_y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_frame
            })

    _windows_dictionary, _windows_to_compute = compute.windows(_arguments,
                                                               _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute, _saturation=True)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple][:2] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _experiments_saturation = {
        _key: [(_fiber_densities[_tuple][2], _fiber_densities[_tuple][1]) for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _tuples_by_experiment = organize.by_experiment(_tuples)

    _correlations_array = []
    _saturation_array = []
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

            # saturation
            _left_cell_saturation = _experiments_saturation[(_experiment, _series, _group, 'left_cell')]
            _right_cell_saturation = _experiments_saturation[(_experiment, _series, _group, 'right_cell')]

            _left_cell_saturation = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['left_cell'],
                _left_cell_saturation
            )
            _right_cell_saturation = compute.remove_blacklist(
                _experiment,
                _series,
                _properties['cells_ids']['right_cell'],
                _right_cell_saturation
            )

            _left_cell_saturation_filtered, _right_cell_saturation_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _left_cell_saturation, _right_cell_saturation
                )

            _left_saturation_fraction_last_time_frame = _left_cell_saturation_filtered[-1]
            _right_saturation_fraction_last_time_frame = _right_cell_saturation_filtered[-1]
            _group_saturation_fraction_mean = \
                (_left_saturation_fraction_last_time_frame + _right_saturation_fraction_last_time_frame) / 2

            _correlations_array.append(_correlation)
            _saturation_array.append(_group_saturation_fraction_mean)

    print('Total points:', len(_correlations_array))
    print('Wilcoxon of correlations around the zero:')
    print(wilcoxon(_correlations_array))
    print('Pearson correlation of correlations and saturation fraction mean:')
    print(compute_lib.correlation(_correlations_array, _saturation_array, _with_p_value=True))

    return _correlations_array, _saturation_array


def main(_offset_y=0.5):
    _correlations_array, _saturation_array = compute_fiber_densities(_offset_y)

    # plot
    _fig = go.Figure(
        data=go.Scatter(
            x=_saturation_array,
            y=_correlations_array,
            mode='markers',
            marker={
                'size': 10,
                'color': '#ea8500'
            },
            showlegend=False
        ),
        layout={
            'xaxis': {
                'title': 'End mean saturation fraction',
                'zeroline': False,
                'range': [-0.01, 0.1],
                'tickmode': 'array',
                'tickvals': [0, 0.05, 0.1]
            },
            'yaxis': {
                'title': 'Inner windows correlation',
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
        _filename='plot_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
