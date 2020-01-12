import random

import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import load, config, filtering, organize, compute
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS

MINIMUM_TIME_POINTS = 20
OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_MIN = 0
CELLS_DISTANCE_MAX = 20
DIRECTION = 'inside'


def main():
    _experiments = load.experiment_groups_as_tuples('SN16_CZI')
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance(
        _experiments, _min_distance=CELLS_DISTANCE_MIN, _max_distance=CELLS_DISTANCE_MAX
    )
    print('num of experiments', len(_experiments))
    random.shuffle(_experiments)
    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in range(len(_experiments)):
        _same_tuple = _experiments[_same_index]
        _same_experiment, _same_series, _same_group = _same_tuple
        _same_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _experiment=_same_experiment,
            _series_id=_same_series,
            _group=_same_group,
            _length_x=ROI_LENGTH,
            _length_y=ROI_WIDTH,
            _length_z=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _offset_z=OFFSET_Z,
            _cell_id='left_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _same_right_cell_fibers_densities = compute.roi_fibers_density_by_time(
            _experiment=_same_experiment,
            _series_id=_same_series,
            _group=_same_group,
            _length_x=ROI_LENGTH,
            _length_y=ROI_WIDTH,
            _length_z=ROI_HEIGHT,
            _offset_x=OFFSET_X,
            _offset_y=OFFSET_Y,
            _offset_z=OFFSET_Z,
            _cell_id='right_cell',
            _direction=DIRECTION,
            _time_points=MINIMUM_TIME_POINTS
        )
        _same_correlation = compute_lib.correlation(
            compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
            compute_lib.derivative(_same_right_cell_fibers_densities, _n=DERIVATIVE)
        )
        for _different_index in range(_same_index + 1, len(_experiments)):
            _different_tuple = _experiments[_different_index]
            _different_experiment, _different_series, _different_group = _different_tuple
            _different_left_cell_fibers_densities = compute.roi_fibers_density_by_time(
                _experiment=_different_experiment,
                _series_id=_different_series,
                _group=_different_group,
                _length_x=ROI_LENGTH,
                _length_y=ROI_WIDTH,
                _length_z=ROI_HEIGHT,
                _offset_x=OFFSET_X,
                _offset_y=OFFSET_Y,
                _offset_z=OFFSET_Z,
                _cell_id='left_cell',
                _direction=DIRECTION,
                _time_points=MINIMUM_TIME_POINTS
            )
            _different_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_same_left_cell_fibers_densities, _n=DERIVATIVE),
                compute_lib.derivative(_different_left_cell_fibers_densities, _n=DERIVATIVE)
            ))
            _same_correlations_array.append(_same_correlation)

    _same_minus_different = np.array(_same_correlations_array) - np.array(_different_correlations_array)
    _same_count = len(_same_minus_different[_same_minus_different > 0])
    _same_percentages = round(_same_count / len(_same_minus_different), 10)
    _wilcoxon_rank = wilcoxon(_same_minus_different)

    print('Same Network:', str(_same_percentages * 100) + '%')
    print('Wilcoxon:', _wilcoxon_rank)


if __name__ == '__main__':
    main()
