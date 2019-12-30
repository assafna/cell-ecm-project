import os

from libs.experiments import paths


def file_data(_file_path):
    try:
        _data = []
        with open(_file_path, 'r') as _file:
            _lines = _file.readlines()
            for _line in _lines[1:]:
                _density_value = float(str(_line.split('\t')[1]).split('\n')[0])
                _data.append(_density_value)
    finally:
        _file.close()

    return _data


def serieses(_experiment_path):
    return [_series for _series in os.listdir(_experiment_path) if
            os.path.isdir(os.path.join(_experiment_path, _series)) and _series.startswith('Series')]


def z_group_fibers_densities(_group_path, _z_group):
    _z_group_path = os.path.join(_group_path, _z_group)
    _z_group_files = [_file_name for _file_name in os.listdir(_z_group_path) if _file_name.startswith('tp_')]
    _tps = {}
    for _file_name in _z_group_files:
        _tp = int(str(_file_name.split('.')[0]).split('_')[1])
        _tp_path = os.path.join(_z_group_path, _file_name)
        _tps[_tp] = file_data(_tp_path)

    return _tps


def group_fibers_densities(_series_path, _group):
    _group_path = os.path.join(_series_path, _group)
    _group_z_groups = [_z_group for _z_group in os.listdir(_group_path) if
                       os.path.isdir(os.path.join(_group_path, _z_group)) and _z_group.startswith('cell')]
    _group_data = {}
    for _z_group in _group_z_groups:
        _group_data[_z_group] = z_group_fibers_densities(_group_path, _z_group)

    return _group_data


def series_fibers_densities(_experiment_path, _series):
    _series_path = os.path.join(_experiment_path, _series)
    _series_groups = [_group for _group in os.listdir(_series_path) if
                      os.path.isdir(os.path.join(_series_path, _group)) and _group.startswith('cell')]
    _series_data = {}
    for _group in _series_groups:
        _series_data[_group] = group_fibers_densities(_series_path, _group)

    return _series_data


def experiment_fibers_densities(_experiment):
    _experiment_path = paths.fibers_density(_experiment)
    _experiment_serieses = serieses(_experiment_path)
    _experiment_data = {}
    for _series in _experiment_serieses:
        _experiment_data[_series] = series_fibers_densities(_experiment_path, _series)

    return _experiment_data


def group_normalization(_series_path, _group):
    _group_path = os.path.join(_series_path, _group)
    return file_data(_group_path)[0]


def series_normalization(_experiment_path, _series):
    _series_path = os.path.join(_experiment_path, _series)
    _series_groups = [_group for _group in os.listdir(_series_path) if _group.startswith('cell')]
    _series_data = {}
    for _group in _series_groups:
        _series_data[str(_group.split('.')[0])] = group_normalization(_series_path, _group)

    return _series_data


def experiment_normalization(_experiment):
    _experiment_path = paths.normalization(_experiment)
    _experiment_serieses = serieses(_experiment_path)
    _experiment_data = {}
    for _series in _experiment_serieses:
        _experiment_data[_series] = series_normalization(_experiment_path, _series)

    return _experiment_data
