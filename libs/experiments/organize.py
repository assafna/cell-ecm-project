from libs.experiments import load, compute


def by_tuples(_experiments_dictionary):
    _tuples = []
    for _experiment in _experiments_dictionary:
        _serieses = _experiments_dictionary[_experiment]
        for _series in _serieses:
            _groups = _serieses[_series]
            for _group in _groups:
                _z_groups = _groups[_group]
                for _z_group in _z_groups:
                    _tuples.append((_experiment, _series, _group, _z_group))

    return _tuples


def by_cells_distance(_experiments_tuples):
    _tuples_by_distance = {}
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _left_cell_coordinates = _group_properties['time_points'][0]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][0]['right_cell']['coordinates']
        _cells_distance = int(round(compute.cells_distance_in_cell_size(
            _experiment=_experiment,
            _series_id=_series_id,
            _cell_1_coordinates=[(_left_cell_coordinates['x'], _left_cell_coordinates['y'], _left_cell_coordinates['z'])],
            _cell_2_coordinates=[(_right_cell_coordinates['x'], _right_cell_coordinates['y'], _right_cell_coordinates['z'])]
        )))
        if _cells_distance in _tuples_by_distance:
            _tuples_by_distance[_cells_distance].append(_tuple)
        else:
            _tuples_by_distance[_cells_distance] = [_tuple]

    return {_distance: _tuples_by_distance[_distance] for _distance in sorted(_tuples_by_distance.keys())}
