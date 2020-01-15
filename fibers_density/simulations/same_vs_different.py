import os
import random

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import compute, filtering, load, paths
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT
from plotting import scatter, save

MINIMUM_TIME_POINTS = 50
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 2
STD = 0.5
CELLS_DISTANCE = 5.0
DIRECTION = 'inside'


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, MINIMUM_TIME_POINTS)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=True,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False
    )
    _simulations = filtering.by_heterogeneity(_simulations, _std=STD)
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in range(len(_simulations)):
        _same_simulation = _simulations[_same_index]
        _same_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _simulation=_same_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='left_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _same_right_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _simulation=_same_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='right_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _same_correlation = compute_lib.correlation(
            compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _different_index in range(_same_index + 1, len(_simulations)):
            _different_simulation = _simulations[_different_index]
            print(_same_simulation, _different_simulation, sep='\t')
            _different_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
                _simulation=_different_simulation,
                _length_x=ROI_WIDTH,
                _length_y=ROI_HEIGHT,
                _offset_x=OFFSET_X,
                _offset_y=OFFSET_Y,
                _cell_id='left_cell',
                _direction=DIRECTION,
                _time_points=MINIMUM_TIME_POINTS
            )
            _different_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_different_left_cell_fibers_densities, _n=DERIVATIVE)
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

    print('Same Network:', str(_same_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create plot


if __name__ == '__main__':
    main()
