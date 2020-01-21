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
    _master_correlations_array = []
    _slave_correlations_array = []
    for _master_index in range(len(_simulations)):
        _master_simulation = _simulations[_master_index]
        _master_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _simulation=_master_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='left_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _master_right_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _simulation=_master_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='right_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _master_correlation = compute_lib.correlation(
            compute_lib.derivative(_master_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_master_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _slave_index in range(_master_index + 1, len(_simulations)):
            _slave_simulation = _simulations[_slave_index]
            print(_master_simulation, _slave_simulation, sep='\t')
            _slave_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
                _simulation=_slave_simulation,
                _length_x=ROI_WIDTH,
                _length_y=ROI_HEIGHT,
                _offset_x=OFFSET_X,
                _offset_y=OFFSET_Y,
                _cell_id='left_cell',
                _direction=DIRECTION,
                _time_points=MINIMUM_TIME_POINTS
            )
            _slave_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_master_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_slave_left_cell_fibers_densities, _n=DERIVATIVE)
            ))
            _master_correlations_array.append(_master_correlation)

    # points plot
    _fig = scatter.create_plot(
        _x_array=[_master_correlations_array],
        _y_array=[_slave_correlations_array],
        _names_array=['Distance ' + str(CELLS_DISTANCE)],
        _modes_array=['markers'],
        _showlegend_array=[False],
        _x_axis_title='master Network Correlation',
        _y_axis_title='slave Network Correlation',
        _title='master vs. slave Network Correlations'
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

    _master_minus_slave = np.array(_master_correlations_array) - np.array(_slave_correlations_array)
    _master_count = len(_master_minus_slave[_master_minus_slave > 0])
    _master_percentages = round(_master_count / len(_master_minus_slave), 10)
    _wilcoxon_rank = wilcoxon(_master_minus_slave)

    print('master Network:', str(_master_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create plot


if __name__ == '__main__':
    main()
