import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import load, filtering, compute, paths
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS
from plotting import scatter, save

MINIMUM_TIME_POINTS = 240
OFFSET_X = (CELL_DIAMETER_IN_MICRONS / 8) * 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 2
CELLS_DISTANCE = 8
DIRECTION = 'inside'


def main():
    _experiments = load.experiment_groups_as_tuples('SN41')
    # _experiments = filtering.by_distance(_experiments, CELLS_DISTANCE)
    # _minimum_time_points = compute.minimum_time_points(_experiments)
    _minimum_time_points = MINIMUM_TIME_POINTS
    _experiments = filtering.by_time_points_amount(_experiments, _minimum_time_points)

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

    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in range(len(_experiments)):
        _same_tuple = _experiments[_same_index]
        _same_experiment, _same_series, _same_group = _same_tuple
        _same_left_cell_fibers_densities = _fibers_densities[(_same_experiment, _same_series, _same_group)]['left_cell']
        _same_right_cell_fibers_densities = _fibers_densities[(_same_experiment, _same_series, _same_group)][
            'right_cell']
        _same_correlation = compute_lib.correlation(
            compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        print(_same_tuple, _same_correlation)
        for _different_index in range(len(_experiments)):
            if _same_index != _different_index:
                _different_tuple = _experiments[_different_index]
                _different_experiment, _different_series, _different_group = _different_tuple
                _different_left_cell_fibers_densities = _fibers_densities[(_different_experiment, _different_series,
                                                                           _different_group)]['left_cell']
                _different_right_cell_fibers_densities = _fibers_densities[(_different_experiment, _different_series,
                                                                           _different_group)]['right_cell']
                _different_correlations_array.append(compute_lib.correlation(
                    compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_different_left_cell_fibers_densities, _n=DERIVATIVE)
                ))
                _same_correlations_array.append(_same_correlation)
                _different_correlations_array.append(compute_lib.correlation(
                    compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_different_right_cell_fibers_densities, _n=DERIVATIVE)
                ))
                _same_correlations_array.append(_same_correlation)
                _different_correlations_array.append(compute_lib.correlation(
                    compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_different_left_cell_fibers_densities, _n=DERIVATIVE)
                ))
                _same_correlations_array.append(_same_correlation)
                _different_correlations_array.append(compute_lib.correlation(
                    compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE),
                    compute_lib.derivative(_different_right_cell_fibers_densities, _n=DERIVATIVE)
                ))
                _same_correlations_array.append(_same_correlation)

    # points plot
    _fig = scatter.create_plot(
        _x_array=[_same_correlations_array],
        _y_array=[_different_correlations_array],
        _names_array=['Distance ' + str(CELLS_DISTANCE)],
        _modes_array=['markers'],
        _showlegend_array=[False],
        _x_axis_title='Same Network Correlation',
        _y_axis_title='Different Network Correlation',
        _title='Same vs. Different Network Correlations'
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
        _filename='plot'
    )

    _same_minus_different = np.array(_same_correlations_array) - np.array(_different_correlations_array)
    _same_count = len(_same_minus_different[_same_minus_different > 0])
    _same_percentages = round(_same_count / len(_same_minus_different), 10)
    _wilcoxon_rank = wilcoxon(_same_minus_different)

    print('N:', str(len(_same_minus_different)))
    print('Same Network:', str(_same_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)


if __name__ == '__main__':
    main()
