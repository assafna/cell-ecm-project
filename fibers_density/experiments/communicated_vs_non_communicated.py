import os
from itertools import product

import numpy as np
from scipy.stats import wilcoxon
from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

EXPERIMENTS = ['SN16']
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 2
CELLS_DISTANCES = [6, 7, 8, 9]
DIRECTION = 'inside'
REAL_CELLS = True
STATIC = False
BAND = True
PLOT = True
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}


def wilcoxon_test(_communicated_correlations_array, _non_communicated_correlations_array, _name):
    _communicated_minus_non_communicated = \
        np.array(_communicated_correlations_array) - np.array(_non_communicated_correlations_array)
    _communicated_count = len(_communicated_minus_non_communicated[_communicated_minus_non_communicated > 0])
    _communicated_percentages = round(_communicated_count / len(_communicated_minus_non_communicated), 10)
    _wilcoxon_rank = wilcoxon(_communicated_minus_non_communicated)

    print(_name)
    print('\tN:', str(len(_communicated_minus_non_communicated)))
    print('\tCommunicated Pair:', str(_communicated_percentages * 100) + '%')
    print('\tWilcoxon:', _wilcoxon_rank)


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': DIRECTION
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {}
    for _key in tqdm(_rois_dictionary):
        _experiments_fibers_densities[_key] = [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]

    _communicated_correlations_array = []
    _non_communicated_correlations_array = []
    for _communicated_index in tqdm(range(len(_experiments)), desc='communicated Loop'):
        _communicated_tuple = _experiments[_communicated_index]
        _communicated_experiment, _communicated_series, _communicated_group = _communicated_tuple

        _communicated_left_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_communicated_experiment, _communicated_series, _communicated_group, 'left_cell')
            ]
        _communicated_right_cell_fibers_densities = \
            _experiments_fibers_densities[
                (_communicated_experiment, _communicated_series, _communicated_group, 'right_cell')
            ]

        _communicated_properties = \
            load.group_properties(_communicated_experiment, _communicated_series, _communicated_group)
        _communicated_left_cell_fibers_densities = compute.remove_blacklist(
            _communicated_experiment,
            _communicated_series,
            _communicated_properties['cells_ids']['left_cell'],
            _communicated_left_cell_fibers_densities
        )
        _communicated_right_cell_fibers_densities = compute.remove_blacklist(
            _communicated_experiment,
            _communicated_series,
            _communicated_properties['cells_ids']['right_cell'],
            _communicated_right_cell_fibers_densities
        )

        _communicated_left_cell_fibers_densities_filtered, _communicated_right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _communicated_left_cell_fibers_densities, _communicated_right_cell_fibers_densities
            )

        # ignore small arrays
        if len(_communicated_left_cell_fibers_densities_filtered) < \
                MINIMUM_CORRELATION_TIME_POINTS[_communicated_experiment]:
            continue

        _communicated_correlation = compute_lib.correlation(
            compute_lib.derivative(_communicated_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_communicated_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
        )
        for _non_communicated_index in range(len(_experiments)):
            if _communicated_index != _non_communicated_index:
                _non_communicated_tuple = _experiments[_non_communicated_index]
                _non_communicated_experiment, _non_communicated_series, _non_communicated_group = \
                    _non_communicated_tuple
                for _communicated_cell_id, _non_communicated_cell_id in product(['left_cell', 'right_cell'],
                                                                                ['left_cell', 'right_cell']):
                    _communicated_fibers_densities = _experiments_fibers_densities[(
                        _communicated_experiment,
                        _communicated_series,
                        _communicated_group,
                        _communicated_cell_id
                    )]
                    _non_communicated_fibers_densities = _experiments_fibers_densities[(
                        _non_communicated_experiment,
                        _non_communicated_series,
                        _non_communicated_group,
                        _non_communicated_cell_id
                    )]

                    _non_communicated_properties = load.group_properties(
                        _non_communicated_experiment, _non_communicated_series, _non_communicated_group
                    )
                    _communicated_fibers_densities = compute.remove_blacklist(
                        _communicated_experiment,
                        _communicated_series,
                        _communicated_properties['cells_ids'][_communicated_cell_id],
                        _communicated_fibers_densities
                    )
                    _non_communicated_fibers_densities = compute.remove_blacklist(
                        _non_communicated_experiment,
                        _non_communicated_series,
                        _non_communicated_properties['cells_ids'][_non_communicated_cell_id],
                        _non_communicated_fibers_densities
                    )

                    _communicated_fibers_densities_filtered, _non_communicated_fibers_densities_filtered = \
                        compute.longest_same_indices_shared_in_borders_sub_array(
                            _communicated_fibers_densities, _non_communicated_fibers_densities
                        )

                    # ignore small arrays
                    if len(_communicated_fibers_densities_filtered) < \
                            MINIMUM_CORRELATION_TIME_POINTS[_non_communicated_experiment]:
                        continue

                    _non_communicated_correlation = compute_lib.correlation(
                        compute_lib.derivative(_communicated_fibers_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_non_communicated_fibers_densities_filtered, _n=DERIVATIVE)
                    )

                    _communicated_correlations_array.append(_communicated_correlation)
                    _non_communicated_correlations_array.append(_non_communicated_correlation)

    # points plot
    if PLOT:
        _fig = scatter.create_plot(
            _x_array=[_communicated_correlations_array],
            _y_array=[_non_communicated_correlations_array],
            _names_array=['Communicated, Non-communicated'],
            _modes_array=['markers'],
            _show_legend_array=[False],
            _x_axis_title='Communicated Pair Correlation',
            _y_axis_title='Non-communicated Pair Correlation',
            _title=None
        )

        _fig = scatter.add_line(
            _fig=_fig,
            _x1=-1, _y1=-1, _x2=1, _y2=1,
            _name='y = x',
            _color='red',
            _showlegend=False
        )

        save.to_html(
            _fig=_fig,
            _path=os.path.join(paths.PLOTS, save.get_module_name()),
            _filename='plot_' + EXPERIMENTS_STR + '_real_' + str(REAL_CELLS) + '_static_' + str(STATIC) + '_band_' +
                      str(BAND)
        )

    # wilcoxon
    if len(_communicated_correlations_array) > 0:
        wilcoxon_test(
            _communicated_correlations_array,
            _non_communicated_correlations_array,
            'Communicated, Non-communicated'
        )


if __name__ == '__main__':
    main()
