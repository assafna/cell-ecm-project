import numpy as np
from scipy.stats import wilcoxon

import computations
import compute_simulation
import filter_simulations
from configurations import ROI_WIDTH, ROI_HEIGHT

TIME_POINTS = 51
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
CELLS_DISTANCE = 3.0


def main():
    _simulations = ['3D_1', '3D_2', '3D_3', '3D_4', '3D_5']
    _simulations = filter_simulations.by_distance(_simulations, CELLS_DISTANCE)
    _insides_correlations_array = []
    _inside_outside_correlations_array = []
    for _simulation in _simulations:
        _left_cell_inside_fibers_densities = compute_simulation.roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='left_cell',
            _direction='inside',
            _time_points=TIME_POINTS
        )
        _left_cell_outside_fibers_densities = compute_simulation.roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='left_cell',
            _direction='outside',
            _time_points=TIME_POINTS
        )
        _right_cell_inside_fibers_densities = compute_simulation.roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='right_cell',
            _direction='inside',
            _time_points=TIME_POINTS
        )
        _right_cell_outside_fibers_densities = compute_simulation.roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=ROI_WIDTH,
            _length_y=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _cell_id='right_cell',
            _direction='outside',
            _time_points=TIME_POINTS
        )
        # left cell
        _insides_correlations_array.append(computations.correlation(
            computations.derivative(_left_cell_inside_fibers_densities, _n=DERIVATIVE),
            computations.derivative(_right_cell_inside_fibers_densities, _n=DERIVATIVE)
        ))
        _inside_outside_correlations_array.append(computations.correlation(
            computations.derivative(_left_cell_inside_fibers_densities, _n=DERIVATIVE),
            computations.derivative(_left_cell_outside_fibers_densities, _n=DERIVATIVE)
        ))
        # right cell
        _insides_correlations_array.append(computations.correlation(
            computations.derivative(_left_cell_inside_fibers_densities, _n=DERIVATIVE),
            computations.derivative(_right_cell_inside_fibers_densities, _n=DERIVATIVE)
        ))
        _inside_outside_correlations_array.append(computations.correlation(
            computations.derivative(_right_cell_inside_fibers_densities, _n=DERIVATIVE),
            computations.derivative(_right_cell_outside_fibers_densities, _n=DERIVATIVE)
        ))

    _insides_minus_inside_outside = np.array(_insides_correlations_array) - np.array(_inside_outside_correlations_array)
    _insides_count = len(_insides_minus_inside_outside[_insides_minus_inside_outside > 0])
    _insides_percentages = round(_insides_count / len(_insides_minus_inside_outside), 10)
    _wilcoxon_rank = wilcoxon(_insides_minus_inside_outside)

    print('Insides:', str(_insides_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create plots


if __name__ == '__main__':
    main()
