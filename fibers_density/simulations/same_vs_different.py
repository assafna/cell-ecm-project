import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import compute, filtering
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT

MINIMUM_TIME_POINTS = 51
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
CELLS_DISTANCE = 3.0
DIRECTION = 'inside'


def main():
    _simulations = ['3D_1', '3D_2', '3D_3', '3D_4', '3D_5']
    _simulations = filtering.by_time_points_amount(_simulations, MINIMUM_TIME_POINTS)
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in range(len(_simulations)):
        for _different_index in range(_same_index + 1, len(_simulations)):
            _same_simulation = _simulations[_same_index]
            _different_simulation = _simulations[_different_index]
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
            _same_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE)
            ))
            _different_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_different_left_cell_fibers_densities, _n=DERIVATIVE)
            ))

    _same_minus_different = np.array(_same_correlations_array) - np.array(_different_correlations_array)
    _same_count = len(_same_minus_different[_same_minus_different > 0])
    _same_percentages = round(_same_count / len(_same_minus_different), 10)
    _wilcoxon_rank = wilcoxon(_same_minus_different)

    print('Same Network:', str(_same_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create plot


if __name__ == '__main__':
    main()
