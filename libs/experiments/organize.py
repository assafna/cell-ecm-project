from libs.experiments import load, compute


def by_pair_distance(_experiments_tuples):
    _tuples_by_distance = {}
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_coordinates = _group_properties['time_points'][0]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][0]['right_cell']['coordinates']
        _pair_distance = int(round(compute.pair_distance_in_cell_size(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_1_coordinates=
            [(_left_cell_coordinates['x'], _left_cell_coordinates['y'], _left_cell_coordinates['z'])],
            _cell_2_coordinates=
            [(_right_cell_coordinates['x'], _right_cell_coordinates['y'], _right_cell_coordinates['z'])]
        )))
        if _pair_distance in _tuples_by_distance:
            _tuples_by_distance[_pair_distance].append(_tuple)
        else:
            _tuples_by_distance[_pair_distance] = [_tuple]

    return {_distance: _tuples_by_distance[_distance] for _distance in sorted(_tuples_by_distance.keys())}


def by_single_cell_id(_experiments_tuples):
    _tuples_by_single_cell_id = {}
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _, _cell_id, _degrees_xy, _degrees_z = _group.split('_')
        _new_tuple = (_experiment, _series_id, _cell_id)
        if _new_tuple in _tuples_by_single_cell_id:
            _tuples_by_single_cell_id[_new_tuple].append(_tuple)
        else:
            _tuples_by_single_cell_id[_new_tuple] = [_tuple]

    return _tuples_by_single_cell_id


def by_matched_real_and_fake(_experiments_tuples):
    _experiments_matched = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _type, _cell_1_id, _cell_2_id = _group.split('_')

        if _type != 'cells':
            continue

        _fake_group = 'fake_' + _cell_1_id + '_' + _cell_2_id
        _fake_tuple = (_experiment, _series_id, _fake_group)
        if _fake_tuple in _experiments_tuples:
            _experiments_matched.append((_tuple, _fake_tuple))

    return _experiments_matched


def by_matched_real_real_and_real_fake(_experiments_tuples, _real_fake_same_cell=True):
    _experiments_matched = []
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_split = _group.split('_')

        if len(_group_split) < 4:
            continue

        _, _cell_1_id, _cell_2_id, _, _cell_real_id, _, _cell_fake_id = _group_split

        if (_cell_real_id == _cell_fake_id) != _real_fake_same_cell:
            continue

        _real_real_group = 'cells_' + _cell_1_id + '_' + _cell_2_id
        _real_real_tuple = (_experiment, _series_id, _real_real_group)

        if _real_real_tuple in _experiments_tuples:
            _experiments_matched.append((_real_real_tuple, _tuple))

    return _experiments_matched


def by_experiment(_experiments_tuples):
    _tuples_by_experiment = {}
    for _tuple in _experiments_tuples:
        _experiment, _, _ = _tuple

        if _experiment in _tuples_by_experiment:
            _tuples_by_experiment[_experiment].append(_tuple)
        else:
            _tuples_by_experiment[_experiment] = [_tuple]

    return _tuples_by_experiment
