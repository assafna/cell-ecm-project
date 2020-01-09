import numpy as np
from scipy.stats import wilcoxon

from libs import compute_lib
from libs.experiments import load, config, filtering, organize, compute

MINIMUM_TIME_POINTS = 15
WINDOW_OFFSET = 5
DERIVATIVE = 1
CELLS_DISTANCE_MIN = 7.0
CELLS_DISTANCE_MAX = 9.0


def main():
    _experiments = load.fibers_density_dictionary(config.PAIRS)
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)
    _experiments = filtering.by_distance(
        _experiments, _min_distance=CELLS_DISTANCE_MIN, _max_distance=CELLS_DISTANCE_MAX
    )
    _experiments = organize.by_tuples(_experiments)
    _same_correlations_array = []
    _different_correlations_array = []
    for _same_index in range(len(_experiments)):
        _same_experiment = _experiments[_same_index]
        _same_experiment_fibers_densities = load.fibers_density_z_group_file_data(
            _same_experiment[0], _same_experiment[1], _same_experiment[2], _same_experiment[3]
        )
        _same_experiment_first_cell = []
        _same_experiment_second_cell = []
        _same_experiment_fibers_densities = compute.fibers_density_cut_edges(_same_experiment_fibers_densities)
        for _time_point in range(MINIMUM_TIME_POINTS):
            _same_experiment_first_cell.append(_same_experiment_fibers_densities[_time_point + 1][WINDOW_OFFSET])
            _same_experiment_second_cell.append(
                _same_experiment_fibers_densities[_time_point + 1][-WINDOW_OFFSET - 1]
            )
        _same_correlation = compute_lib.correlation(
            compute_lib.derivative(_same_experiment_first_cell, _n=DERIVATIVE),
            compute_lib.derivative(_same_experiment_second_cell, _n=DERIVATIVE)
        )
        for _different_index in range(_same_index + 1, len(_experiments)):
            _different_experiment = _experiments[_different_index]
            _different_experiment_fibers_densities = load.fibers_density_z_group_file_data(
                _different_experiment[0], _different_experiment[1], _different_experiment[2], _different_experiment[3]
            )
            _different_experiment_first_cell = []
            _different_experiment_fibers_densities = compute.fibers_density_cut_edges(
                _different_experiment_fibers_densities
            )
            for _time_point in range(MINIMUM_TIME_POINTS):
                _different_experiment_first_cell.append(
                    _different_experiment_fibers_densities[_time_point + 1][WINDOW_OFFSET]
                )
            _different_correlations_array.append(compute_lib.correlation(
                compute_lib.derivative(_same_experiment_first_cell, _n=DERIVATIVE),
                compute_lib.derivative(_different_experiment_first_cell, _n=DERIVATIVE)
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
