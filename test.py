from libs.compute_lib import distance_from_a_point_to_a_line
from libs.simulations import load, filtering, compute


if __name__ == '__main__':
    # _experiments = load.experiment_groups_as_tuples(_experiment='SN16')
    # _experiments = filtering.by_time_points_amount(_experiments, _time_points=18)
    # _experiments = filtering.by_distances(_experiments, _distances=[6, 7, 8, 9])
    # _experiments = filtering.by_real_cells(_experiments)
    # _experiments = filtering.by_band(_experiments)
    # print(len(_experiments))
    # _index = 0
    # print('Index', 'Start distance', 'End distance', 'Difference', sep='\t')
    # for _tuple in _experiments:
    #     _experiment, _series_id, _group = _tuple
    #     _distance_start = compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point=0)
    #     _distance_end = compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point=17)
    #     print(_index, _distance_start, _distance_end, _distance_start - _distance_end, sep='\t')
    #     _index += 1
    # _a = distance_from_a_point_to_a_line([-1, -1, 1, 1], [1, 4])
    # print(_a)
    _simulations = load.structured()
    # _simulations = filtering.by_categories(
    #     _simulations,
    #     _is_single_cell=False,
    #     _is_heterogeneity=True,
    #     _is_low_connectivity=False,
    #     _is_causality=False,
    #     _is_dominant_passive=True
    # )
    _properties = load.properties(_simulation='7D_8_heterogeneity')
    print(_properties)
    # print(_simulations)

