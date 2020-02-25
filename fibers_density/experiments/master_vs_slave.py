import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT
from methods.experiments import export_video
from plotting import scatter, save

EXPERIMENTS = ['SN16']
EXPERIMENTS_STR = '_'.join(EXPERIMENTS)
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 2
CELLS_DISTANCES = [6, 7, 8]
DIRECTION = 'inside'
REAL_CELLS = True
STATIC = False
BAND = True
PLOT = False
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50
}


def wilcoxon_test(_master_correlations_array, _slave_correlations_array, _name):
    _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
    _master_count = len(_master_minus_slave[_master_minus_slave > 0])
    _master_percentages = round(_master_count / len(_master_minus_slave), 10)
    _wilcoxon_rank = wilcoxon(_master_minus_slave)

    print(_name)
    print('\tN:', str(len(_master_minus_slave)))
    print('\tMaster Network:', str(_master_percentages * 100) + '%')
    print('\tWilcoxon:', _wilcoxon_rank)


def main():
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS)
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    if BAND:
        _experiments = filtering.by_band(_experiments)

    print(len(_experiments))

    # prepare data in mp
    _fibers_densities = {}
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series, _group = _tuple
        _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH,
                           OFFSET_X, OFFSET_Y, OFFSET_Z, DIRECTION))

    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
    _p.close()
    for _index in range(len(_arguments)):
        _tuple = _experiments[_index]
        _answer = _answers[_index]
        _fibers_densities[_tuple] = _answer

    _master_correlations_array = []
    _slave_correlations_array = []
    for _master_index in range(len(_experiments)):
        _master_tuple = _experiments[_master_index]
        _master_experiment, _master_series, _master_group = _master_tuple
        _master_left_cell_fibers_densities = \
            _fibers_densities[(_master_experiment, _master_series, _master_group)]['left_cell']
        _master_right_cell_fibers_densities = \
            _fibers_densities[(_master_experiment, _master_series, _master_group)]['right_cell']

        _master_properties = load.group_properties(_master_experiment, _master_series, _master_group)
        _master_left_cell_fibers_densities = compute.remove_blacklist(
            _master_experiment, _master_series, _master_properties['cells_ids']['left_cell'], _master_left_cell_fibers_densities)
        _master_right_cell_fibers_densities = compute.remove_blacklist(
            _master_experiment, _master_series, _master_properties['cells_ids']['right_cell'], _master_right_cell_fibers_densities)

        _master_left_cell_fibers_densities_filtered, _master_right_cell_fibers_densities_filtered = \
            compute.longest_same_indices_shared_in_borders_sub_array(
                _master_left_cell_fibers_densities, _master_right_cell_fibers_densities
            )

        # ignore small arrays
        if len(_master_left_cell_fibers_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_master_experiment]:
            continue

        _master_correlation = compute_lib.correlation(
            compute_lib.derivative(_master_left_cell_fibers_densities_filtered, _n=DERIVATIVE),
            compute_lib.derivative(_master_right_cell_fibers_densities_filtered, _n=DERIVATIVE)
        )
        for _slave_index in range(len(_experiments)):
            if _master_index != _slave_index:
                _slave_tuple = _experiments[_slave_index]
                _slave_experiment, _slave_series, _slave_group = _slave_tuple
                for _master_cell_id, _slave_cell_id in product(['left_cell', 'right_cell'],
                                                               ['left_cell', 'right_cell']):
                    _master_fibers_densities = \
                        _fibers_densities[(_master_experiment, _master_series, _master_group)][_master_cell_id]
                    _slave_fibers_densities = \
                        _fibers_densities[(_slave_experiment, _slave_series, _slave_group)][_slave_cell_id]

                    _slave_properties = load.group_properties(_slave_experiment, _slave_series, _slave_group)
                    _master_fibers_densities = compute.remove_blacklist(
                        _master_experiment, _master_series, _master_properties['cells_ids'][_master_cell_id], _master_fibers_densities)
                    _slave_fibers_densities = compute.remove_blacklist(
                        _slave_experiment, _slave_series, _slave_properties['cells_ids'][_slave_cell_id], _slave_fibers_densities)

                    _master_fibers_densities_filtered, _slave_fibers_densities_filtered = \
                        compute.longest_same_indices_shared_in_borders_sub_array(
                            _master_fibers_densities, _slave_fibers_densities
                        )

                    # ignore small arrays
                    if len(_master_fibers_densities_filtered) < MINIMUM_CORRELATION_TIME_POINTS[_slave_experiment]:
                        continue

                    _slave_correlation = compute_lib.correlation(
                        compute_lib.derivative(_master_fibers_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_slave_fibers_densities_filtered, _n=DERIVATIVE)
                    )

                    _master_correlations_array.append(_master_correlation)
                    _slave_correlations_array.append(_slave_correlation)

    # points plot
    if PLOT:
        _fig = scatter.create_plot(
            _x_array=[_master_correlations_array],
            _y_array=[_slave_correlations_array],
            _names_array=['Master, Slave'],
            _modes_array=['markers'],
            _showlegend_array=[False],
            _x_axis_title='Master Network Correlation',
            _y_axis_title='Slave Network Correlation',
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
    if len(_master_correlations_array) > 0:
        wilcoxon_test(_master_correlations_array, _slave_correlations_array, 'Master, Slave')


if __name__ == '__main__':
    main()
