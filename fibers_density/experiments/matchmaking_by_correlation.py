from tqdm import tqdm

from libs import compute_lib
from libs.experiments import load, filtering, compute
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT

# based on time resolution
EXPERIMENTS = {
    False: ['SN16'],
    True: ['SN41']
}
OFFSET_X = 0
# TODO: set the offset in y according to the angle in the original Z slices of the cells
OFFSET_Y = 0.5
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCE_RANGE = [4, 10]
REAL_CELLS = True
STATIC = False
MINIMUM_CORRELATION_TIME_POINTS = {
    'SN16': 15,
    'SN18': 15,
    'SN41': 50,
    'SN44': 50,
    'SN45': 50
}


def main(_band=True, _high_time_resolution=False):
    _experiments = load.experiments_groups_as_tuples(EXPERIMENTS[_high_time_resolution])
    _experiments = filtering.by_distance_range(_experiments, CELLS_DISTANCE_RANGE)
    _experiments = filtering.by_real_cells(_experiments, _real_cells=REAL_CELLS)
    _experiments = filtering.by_static_cells(_experiments, _static=STATIC)
    _experiments = filtering.by_band(_experiments, _band=_band)
    print('Total experiments:', len(_experiments))

    _arguments = []
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple

        # stop when windows are overlapping
        _properties = load.group_properties(_experiment, _series_id, _group)
        _latest_time_point = len(_properties['time_points'])
        for _time_point in range(len(_properties['time_points'])):
            _cells_distance = \
                compute.cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point)
            if _cells_distance - 1 - OFFSET_X * 2 < ROI_LENGTH * 2:
                _latest_time_point = _time_point - 1
                break

        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'experiment': _experiment,
                'series_id': _series_id,
                'group': _group,
                'length_x': ROI_LENGTH,
                'length_y': ROI_HEIGHT,
                'length_z': ROI_WIDTH,
                'offset_x': OFFSET_X,
                'offset_y': OFFSET_Y,
                'offset_z': OFFSET_Z,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': _latest_time_point
            })

    _rois_dictionary, _rois_to_compute = compute.rois(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id'])
    _fibers_densities = compute.fibers_densities(_rois_to_compute)

    _experiments_fibers_densities = {
        _key: [_fibers_densities[_tuple] for _tuple in _rois_dictionary[_key]]
        for _key in _rois_dictionary
    }

    _correct_matches = 0
    _incorrect_matches = 0
    _n_total = 0
    for _tuple_1 in tqdm(_experiments, desc='Main loop'):
        _experiment_1, _series_id_1, _group_1 = _tuple_1
        _experiment_1_properties = load.group_properties(_experiment_1, _series_id_1, _group_1)

        for _cell_1_id, _real_matched_cell_id in zip(['left_cell', 'right_cell'], ['right_cell', 'left_cell']):
            _cell_1_fibers_densities = \
                _experiments_fibers_densities[(_experiment_1, _series_id_1, _group_1, _cell_1_id)]
            _cell_1_fibers_densities = compute.remove_blacklist(
                _experiment_1,
                _series_id_1,
                _experiment_1_properties['cells_ids'][_cell_1_id],
                _cell_1_fibers_densities
            )

            _best_match = None
            _highest_correlation = -1
            _n = 0
            for _tuple_2 in _experiments:
                _experiment_2, _series_id_2, _group_2 = _tuple_2
                _experiment_2_properties = load.group_properties(_experiment_2, _series_id_2, _group_2)

                for _cell_2_id in ['left_cell', 'right_cell']:

                    # same cell
                    if _tuple_1 == _tuple_2 and _cell_1_id == _cell_2_id:
                        continue

                    _cell_2_fibers_densities = \
                        _experiments_fibers_densities[(_experiment_2, _series_id_2, _group_2, _cell_2_id)]
                    _cell_2_fibers_densities = compute.remove_blacklist(
                        _experiment_2,
                        _series_id_2,
                        _experiment_2_properties['cells_ids'][_cell_2_id],
                        _cell_2_fibers_densities
                    )

                    _cell_1_fibers_densities_filtered, _cell_2_fibers_densities_filtered = \
                        compute.longest_same_indices_shared_in_borders_sub_array(
                            _cell_1_fibers_densities, _cell_2_fibers_densities
                        )

                    # ignore small arrays
                    if len(_cell_1_fibers_densities_filtered) < \
                            MINIMUM_CORRELATION_TIME_POINTS[_experiment_1]:
                        continue

                    _correlation = compute_lib.correlation(
                        compute_lib.derivative(_cell_1_fibers_densities_filtered, _n=DERIVATIVE),
                        compute_lib.derivative(_cell_2_fibers_densities_filtered, _n=DERIVATIVE)
                    )

                    if _correlation > _highest_correlation:
                        _highest_correlation = _correlation
                        _best_match = (_tuple_2, _cell_2_id)

                    _n += 1

            # check matchmaking
            if _best_match is not None:
                if _best_match == (_tuple_1, _real_matched_cell_id):
                    _correct_matches += 1
                else:
                    _incorrect_matches += 1
                _n_total += _n

    _total_cells = _correct_matches + _incorrect_matches
    print('Matchmaking results:')
    print('Total cells:', _total_cells)
    print('Average potential matches per cell:', round(_n_total / _total_cells, 2))
    print('Correct match probability:', round(1 / (_n_total / _total_cells), 2))
    print('Fraction of correct matches:', round(_correct_matches / _total_cells, 2))


if __name__ == '__main__':
    main()
