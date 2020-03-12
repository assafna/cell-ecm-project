import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.simulations import compute, filtering
from libs.simulations.config import ROI_WIDTH, ROI_HEIGHT

TIME_POINTS = 51
OFFSET_X = 0
OFFSET_Y = 0
DERIVATIVE = 1
CELLS_DISTANCE = 3.0


def main():
    _simulations = ['3D_1', '3D_2', '3D_3', '3D_4', '3D_5']
    _simulations = filtering.by_distance(_simulations, _distance=CELLS_DISTANCE)
    _insides_correlations_array = []
    _outsides_correlations_array = []
    for _simulation in _simulations:
        _left_cell_inside_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'left_cell',
            'direction': 'inside',
            'time_points': TIME_POINTS
        })
        _left_cell_outside_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'left_cell',
            'direction': 'outside',
            'time_points': TIME_POINTS
        })
        _right_cell_inside_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'right_cell',
            'direction': 'inside',
            'time_points': TIME_POINTS
        })
        _right_cell_outside_fibers_densities = compute.roi_fibers_density_by_time({
            'simulation': _simulation,
            'length_x': ROI_WIDTH,
            'length_y': ROI_HEIGHT,
            'offset_x': OFFSET_X,
            'offset_y': OFFSET_Y,
            'cell_id': 'right_cell',
            'direction': 'outside',
            'time_points': TIME_POINTS
        })
        _insides_correlations_array.append(compute_lib.correlation(
            compute_lib.derivative(_left_cell_inside_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_inside_fibers_densities, _n=DERIVATIVE)
        ))
        _outsides_correlations_array.append(compute_lib.correlation(
            compute_lib.derivative(_left_cell_outside_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_right_cell_outside_fibers_densities, _n=DERIVATIVE)
        ))

    _insides_minus_outsides = np.array(_insides_correlations_array) - np.array(_outsides_correlations_array)
    _insides_count = len(_insides_minus_outsides[_insides_minus_outsides > 0])
    _insides_percentages = round(_insides_count / len(_insides_minus_outsides), 10)
    _wilcoxon_rank = wilcoxon(_insides_minus_outsides)

    print('Insides:', str(_insides_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)

    # TODO: create plots


if __name__ == '__main__':
    main()
