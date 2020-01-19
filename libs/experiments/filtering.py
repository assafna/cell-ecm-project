from libs.experiments import load, compute, paths


def by_distance(_experiments_tuples, _distance):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _cell_coordinates_tracked_file = 'series_' + str(_series_id) + '.txt'
        _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(
            _experiment, _cell_coordinates_tracked_file
        )
        _cell_1_id, _cell_2_id = _group.split('_')[1], _group.split('_')[2]
        _cell_1_coordinates = _cell_coordinates_tracked[int(_cell_1_id)]
        _cell_2_coordinates = _cell_coordinates_tracked[int(_cell_2_id)]
        _cells_distance = compute.cells_distance_in_cell_size(
            _experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates
        )
        if _distance - 0.5 <= _cells_distance <= _distance + 0.5:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_distances(_experiments_tuples, _distances):
    # TODO: improve efficiency
    _experiments_tuples_filtered = []
    for _distance in _distances:
        _distance_tuples = by_distance(_experiments_tuples, _distance)
        for _tuple in _distance_tuples:
            if _tuple not in _experiments_tuples_filtered:
                _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_time_points_amount(_experiments_tuples, _time_points, _exactly=False):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _time_points_amount = len(_group_properties['time_points'])
        if (_exactly and _time_points_amount == _time_points) or \
                (not _exactly and _time_points_amount >= _time_points):
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_band(_experiments_tuples, _band=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _experiments_tuples_filtered.append(_tuple) if _group_properties['band'] == _band else None

    return _experiments_tuples_filtered
