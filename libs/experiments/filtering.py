from libs.experiments import load, compute, paths


def by_distance(_experiments_tuples, _distance, _time_point=0):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if len(_group_properties['time_points']) <= _time_point:
            continue
        _left_cell_coordinates = _group_properties['time_points'][_time_point]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][_time_point]['right_cell']['coordinates']
        _cells_distance = compute.cells_distance_in_cell_size(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_1_coordinates=[(_left_cell_coordinates['x'], _left_cell_coordinates['y'], _left_cell_coordinates['z'])],
            _cell_2_coordinates=[(_right_cell_coordinates['x'], _right_cell_coordinates['y'], _right_cell_coordinates['z'])]
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


def by_series_id(_experiments_tuples, _series_id):
    return [_tuple for _tuple in _experiments_tuples if _tuple[1] == _series_id]


def by_real_cells(_experiments_tuples, _real_cells=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if _group_properties['fake'] != _real_cells:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered


def by_static_cells(_experiments_tuples, _static=True):
    _experiments_tuples_filtered = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        if _group_properties['static'] == _static:
            _experiments_tuples_filtered.append(_tuple)

    return _experiments_tuples_filtered
