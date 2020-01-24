import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS
from plotting import scatter, save

MINIMUM_TIME_POINTS = 15
OFFSET_X = CELL_DIAMETER_IN_MICRONS * 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = CELL_DIAMETER_IN_MICRONS * 1
OFFSET_Z = CELL_DIAMETER_IN_MICRONS * 0
DERIVATIVE = 2
CELLS_DISTANCES = [8]
DIRECTION = 'inside'


def wilcoxon_test(_master_correlations_array, _slave_correlations_array, _name):
    _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
    _master_count = len(_master_minus_slave[_master_minus_slave > 0])
    _master_percentages = round(_master_count / len(_master_minus_slave), 10)
    _wilcoxon_rank = wilcoxon(_master_minus_slave)

    print(_name)
    print('\tN:', str(len(_master_minus_slave)))
    print('\tmaster Network:', str(_master_percentages * 100) + '%')
    print('\tWilcoxon:', _wilcoxon_rank)


def main():
    _experiments = load.experiment_groups_as_tuples('SN16')
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments_band = filtering.by_band(_experiments, _band=True)
    # _experiments = _experiments_band
    _minimum_time_points = MINIMUM_TIME_POINTS
    _experiments = filtering.by_time_points_amount(_experiments, _minimum_time_points)
    _experiments.remove(('SN16', 13, 'cells_0_1'))
    # _experiments.remove(('SN16', 21, 'cells_0_1'))
    # _experiments.remove(('SN16', 1, 'cells_1_5'))
    # _experiments = filtering.by_series_id(_experiments, _series_id=22)
    # _experiments = [
    #     ('SN16', 1, 'cells_2_3'),
    #     ('SN16', 1, 'cells_2_4'),
    #     ('SN16', 1, 'cells_2_5'),
    #     # ('SN16', 22, 'cells_0_1')
    # ]

    # prepare data in mp
    _fibers_densities = {}
    _arguments = []
    for _tuple in _experiments:
        _experiment, _series, _group = _tuple
        _arguments.append((_experiment, _series, _group, ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT,
                           OFFSET_X, OFFSET_Y, OFFSET_Z, DIRECTION, _minimum_time_points))
    _p = Pool()
    _answers = _p.starmap(compute.roi_fibers_density_by_time_pairs, _arguments)
    _p.close()
    for _index in range(len(_arguments)):
        _tuple = _experiments[_index]
        _answer = _answers[_index]
        _fibers_densities[_tuple] = _answer

    _band_band_master_correlations_array = []
    _band_band_slave_correlations_array = []
    _band_no_band_master_correlations_array = []
    _band_no_band_slave_correlations_array = []
    _no_band_band_master_correlations_array = []
    _no_band_band_slave_correlations_array = []
    _no_band_no_band_master_correlations_array = []
    _no_band_no_band_slave_correlations_array = []

    for _master_index in range(len(_experiments)):
        _master_tuple = _experiments[_master_index]
        _master_experiment, _master_series, _master_group = _master_tuple
        _master_left_cell_fibers_densities = _fibers_densities[(_master_experiment, _master_series, _master_group)]['left_cell']
        _master_right_cell_fibers_densities = _fibers_densities[(_master_experiment, _master_series, _master_group)][
            'right_cell']
        _master_correlation = compute_lib.correlation(
            compute_lib.derivative(_master_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_master_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        print(_master_tuple, _master_correlation)
        for _slave_index in range(len(_experiments)):
            if _master_index != _slave_index:
                _slave_tuple = _experiments[_slave_index]
                _slave_experiment, _slave_series, _slave_group = _slave_tuple
                _slave_left_cell_fibers_densities = _fibers_densities[(_slave_experiment, _slave_series,
                                                                           _slave_group)]['left_cell']
                _slave_right_cell_fibers_densities = _fibers_densities[(_slave_experiment, _slave_series,
                                                                           _slave_group)]['right_cell']
                _slave_correlation = compute_lib.correlation(
                    compute_lib.derivative(_master_left_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_slave_left_cell_fibers_densities, _n=DERIVATIVE)
                )

                if _master_tuple in _experiments_band and _slave_tuple in _experiments_band:
                    _band_band_master_correlations_array.append(_master_correlation)
                    _band_band_slave_correlations_array.append(_slave_correlation)
                elif _master_tuple in _experiments_band and _slave_tuple not in _experiments_band:
                    _band_no_band_master_correlations_array.append(_master_correlation)
                    _band_no_band_slave_correlations_array.append(_slave_correlation)
                elif _master_tuple not in _experiments_band and _slave_tuple in _experiments_band:
                    _no_band_band_master_correlations_array.append(_master_correlation)
                    _no_band_band_slave_correlations_array.append(_slave_correlation)
                else:
                    _no_band_no_band_master_correlations_array.append(_master_correlation)
                    _no_band_no_band_slave_correlations_array.append(_slave_correlation)

    # points plot
    _fig = scatter.create_plot(
        _x_array=[
            _band_band_master_correlations_array,
            _band_no_band_master_correlations_array,
            _no_band_band_master_correlations_array,
            _no_band_no_band_master_correlations_array
        ],
        _y_array=[
            _band_band_slave_correlations_array,
            _band_no_band_slave_correlations_array,
            _no_band_band_slave_correlations_array,
            _no_band_no_band_slave_correlations_array
        ],
        _names_array=[
            'Master Band, Slave Band',
            'Master Band, Slave No Band',
            'Master No Band, Slave Band',
            'Master No Band, Slave No Band'
        ],
        _modes_array=['markers'] * 4,
        _showlegend_array=[True] * 4,
        _x_axis_title='Master Network Correlation',
        _y_axis_title='Slave Network Correlation',
        _title='Master vs. Slave Network Correlations'
    )

    _fig = scatter.add_line(
        _fig=_fig,
        _x1=-1, _y1=-1, _x2=1, _y2=1,
        _name='y = x',
        _color='red',
        _showlegend=True
    )

    save.to_html(
        _fig=_fig,
        _path=os.path.join(paths.PLOTS, save.get_module_name()),
        _filename='plot_derivative_' + str(DERIVATIVE)
    )

    # wilcoxon
    wilcoxon_test(_band_band_master_correlations_array, _band_band_slave_correlations_array, 'Master Band, Slave Band')
    wilcoxon_test(_band_no_band_master_correlations_array, _band_no_band_slave_correlations_array, 'Master Band, Slave No Band')
    wilcoxon_test(_no_band_band_master_correlations_array, _no_band_band_slave_correlations_array, 'Master No Band, Slave Band')
    wilcoxon_test(_no_band_no_band_master_correlations_array, _no_band_no_band_slave_correlations_array, 'Master No Band, Slave No Band')


if __name__ == '__main__':
    main()
