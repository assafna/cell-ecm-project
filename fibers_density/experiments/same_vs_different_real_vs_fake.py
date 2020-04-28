import os
from itertools import product

import plotly.graph_objs as go
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, organize, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT

# based on time resolution
from plotting import save

EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41']
}
OFFSET_X = 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}


def compute_fibers_densities(_offset_y=0.5, _high_time_resolution=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_band(_experiments)
    _experiments_matched = organize.by_matched_real_and_fake(_experiments)
    print('Total matched pairs:', len(_experiments_matched))

    _arguments = []
    for _matched_tuple in _experiments_matched:
        for _tuple in _matched_tuple:
            _experiment, _series_id, _group = _tuple

            # stop when windows are overlapping
            _properties = load.group_properties(_experiment, _series_id, _group)
            _latest_time_point = len(_properties['time_points'])
            for _time_point in range(len(_properties['time_points'])):
                _cells_distance = \
                    compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
                if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                    _latest_time_point = _time_point - 1
                    break

            for _cell_id in ['left_cell', 'right_cell']:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': ROI_LENGTH,
                    'length_y': ROI_HEIGHT,
                    'length_z': ROI_WIDTH,
                    'offset_x': OFFSET_X,
                    'offset_y': _offset_y,
                    'offset_z': OFFSET_Z,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_points': _latest_time_point
                })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    # same (real, fake), different (real, fake)
    _correlations = [
        [
            [], []
        ],
        [
            [], []
        ]
    ]
    for _same_index in tqdm(range(len(_experiments_matched)), desc='Main loop'):
        for _group_type_index in [0, 1]:
            _same_tuple = _experiments_matched[_same_index][_group_type_index]
            _same_experiment, _same_series, _same_group = _same_tuple

            _same_left_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_same_experiment, _same_series, _same_group, 'left_cell')
                ]
            _same_right_cell_fibers_densities = \
                _experiments_fibers_densities[
                    (_same_experiment, _same_series, _same_group, 'right_cell')
                ]

            _same_properties = \
                load.group_properties(_same_experiment, _same_series, _same_group)
            _same_left_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['left_cell'],
                _same_left_cell_fibers_densities
            )
            _same_right_cell_fibers_densities = compute.remove_blacklist(
                _same_experiment,
                _same_series,
                _same_properties['cells_ids']['right_cell'],
                _same_right_cell_fibers_densities
            )

            _same_left_cell_fibers_densities_filtered, _same_right_cell_fibers_densities_filtered = \
                compute.longest_same_indices_shared_in_borders_sub_array(
                    _same_left_cell_fibers_densities, _same_right_cell_fibers_densities
                )

            # ignore small arrays
            if len(_same_left_cell_fibers_densities_filtered) < \
                    MINIMUM_CORRELATION_TIME_POINTS[_same_experiment]:
                for _different_index in range(len(_experiments_matched)):
                    if _same_index != _different_index:
                        # for all combinations
                        for _i in range(4):
                            _correlations[0][_group_type_index].append(None)
                            _correlations[1][_group_type_index].append(None)
                continue

            _same_correlation = compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
            )
            for _different_index in range(len(_experiments_matched)):
                if _same_index != _different_index:
                    _different_tuple = _experiments_matched[_different_index][_group_type_index]
                    _different_experiment, _different_series, _different_group = _different_tuple
                    for _same_cell_id, _different_cell_id in product(['left_cell', 'right_cell'],
                                                                     ['left_cell', 'right_cell']):
                        _same_fibers_densities = _experiments_fibers_densities[(
                            _same_experiment,
                            _same_series,
                            _same_group,
                            _same_cell_id
                        )]
                        _different_fibers_densities = _experiments_fibers_densities[(
                            _different_experiment,
                            _different_series,
                            _different_group,
                            _different_cell_id
                        )]

                        _different_properties = load.group_properties(
                            _different_experiment, _different_series, _different_group
                        )
                        _same_fibers_densities = compute.remove_blacklist(
                            _same_experiment,
                            _same_series,
                            _same_properties['cells_ids'][_same_cell_id],
                            _same_fibers_densities
                        )
                        _different_fibers_densities = compute.remove_blacklist(
                            _different_experiment,
                            _different_series,
                            _different_properties['cells_ids'][_different_cell_id],
                            _different_fibers_densities
                        )

                        _same_fibers_densities_filtered, _different_fibers_densities_filtered = \
                            compute.longest_same_indices_shared_in_borders_sub_array(
                                _same_fibers_densities, _different_fibers_densities
                            )

                        # ignore small arrays
                        if len(_same_fibers_densities_filtered) < \
                                MINIMUM_CORRELATION_TIME_POINTS[_different_experiment]:
                            _correlations[0][_group_type_index].append(None)
                            _correlations[1][_group_type_index].append(None)
                            continue

                        _different_correlation = compute_lib.correlation(
                            compute_lib.derivative(_same_fibers_densities_filtered, _n=DERIVATIVE),
                            compute_lib.derivative(_different_fibers_densities_filtered, _n=DERIVATIVE)
                        )

                        _correlations[0][_group_type_index].append(_same_correlation)
                        _correlations[1][_group_type_index].append(_different_correlation)

    _distances_from_y_equal_x = [[], []]
    _same_correlations, _different_correlations = _correlations
    _same_real_correlations, _same_fake_correlations = _same_correlations
    _different_real_correlations, _different_fake_correlations = _different_correlations
    for _same_real, _same_fake, _different_real, _different_fake in \
            zip(_same_real_correlations, _same_fake_correlations,
                _different_real_correlations, _different_fake_correlations):

        # one of the correlations is none - not valid
        if None in [_same_real, _same_fake, _different_real, _different_fake]:
            continue

        for _group_type_index, _same, _different in \
                zip([0, 1], [_same_real, _same_fake], [_different_real, _different_fake]):

            _point_distance = compute_lib.distance_from_a_point_to_a_line(
                _line=[-1, -1, 1, 1],
                _point=[_same, _different]
            )
            if _same > _different:
                _distances_from_y_equal_x[_group_type_index].append(_point_distance)
            else:
                _distances_from_y_equal_x[_group_type_index].append(-_point_distance)

    return _distances_from_y_equal_x


def main(_offset_y=0.5, _high_time_resolution=False):
    _distances_from_y_equal_x = compute_fibers_densities(_offset_y, _high_time_resolution)

    print('Total points:', len(_distances_from_y_equal_x[0]))
    print('Wilcoxon between real points and fake points distances from y = x:')
    print(wilcoxon(_distances_from_y_equal_x[0], _distances_from_y_equal_x[1]))

    # plot
    _colors_array = ['#844b00', '#edbc80']
    _fig = go.Figure(
        data=[
            go.Box(
                y=_y_array,
                name=_real_pairs,
                boxpoints=False,
                line={
                    'width': 1
                },
                showlegend=False
            ) for _y_array, _real_pairs, _color in zip(_distances_from_y_equal_x, [True, False], _colors_array)
        ],
        layout={
            'xaxis': {
                'title': 'Real pairs',
                'zeroline': False
            },
            'yaxis': {
                'title': 'Distance from y = x',
                'zeroline': False,
                'range': [-1, 1],
                'tickmode': 'array',
                'tickvals': [-1, -0.5, 0, 0.5, 1]
            }
        }
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_high_time_res_' + str(_high_time_resolution) + '_offset_y_' + str(_offset_y)
    )


if __name__ == '__main__':
    main()
