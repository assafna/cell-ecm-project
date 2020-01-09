from libs.experiments import load, compute, paths


def by_distance(_experiments_dictionary, _min_distance, _max_distance):
    _experiments_dictionary_filtered = {}
    for _experiment in _experiments_dictionary:
        _serieses = _experiments_dictionary[_experiment]
        for _series in _serieses:
            _series_id = str(_series.split()[1])
            _groups = _serieses[_series]
            for _group in _groups:
                _group_id = _group.split('cells_')[1]
                _cell_coordinates_tracked_file = 'series_' + _series_id + '.txt'
                _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(
                    _experiment, _cell_coordinates_tracked_file
                )
                _cell_1_id, _cell_2_id = _group_id.split('_')
                _cell_1_coordinates = _cell_coordinates_tracked[int(_cell_1_id) - 1]
                _cell_2_coordinates = _cell_coordinates_tracked[int(_cell_2_id) - 1]
                _distance = compute.cells_distance_in_cell_size(
                    _experiment, _series, _cell_1_coordinates, _cell_2_coordinates
                )

                if _min_distance <= _distance <= _max_distance:
                    if _experiment not in _experiments_dictionary_filtered:
                        _experiments_dictionary_filtered[_experiment] = {}
                    if _series not in _experiments_dictionary_filtered[_experiment]:
                        _experiments_dictionary_filtered[_experiment][_series] = {}
                    _experiments_dictionary_filtered[_experiment][_series][_group] = _groups[_group]

    return _experiments_dictionary_filtered


def by_time_points_amount(_experiments_dictionary, _time_points, _exactly=False):
    _experiments_dictionary_filtered = {}
    for _experiment in _experiments_dictionary:
        _serieses = _experiments_dictionary[_experiment]
        for _series in _serieses:
            _series_id = str(_series.split()[1])
            _groups = _serieses[_series]
            for _group in _groups:
                _z_groups = _groups[_group]
                for _z_group in _z_groups:
                    _time_points_amount = len(
                        paths.text_files(paths.fibers_density(_experiment, _series, _group, _z_group))
                    )
                    if (_exactly and _time_points_amount == _time_points) or\
                            (not _exactly and _time_points_amount >= _time_points):
                        if _experiment not in _experiments_dictionary_filtered:
                            _experiments_dictionary_filtered[_experiment] = {}
                        if _series not in _experiments_dictionary_filtered[_experiment]:
                            _experiments_dictionary_filtered[_experiment][_series] = {}
                        if _group not in _experiments_dictionary_filtered[_experiment][_series]:
                            _experiments_dictionary_filtered[_experiment][_series][_group] = []
                        _experiments_dictionary_filtered[_experiment][_series][_group].append(_z_group)

    return _experiments_dictionary_filtered
