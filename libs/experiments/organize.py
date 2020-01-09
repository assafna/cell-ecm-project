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
