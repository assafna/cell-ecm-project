import os

from libs.experiments import paths


def objects_time_point_file_data(_experiment, _series, _time_point):
    _file_path = paths.objects(_experiment, _series, _time_point)
    try:
        _cells = []
        with open(_file_path) as _f:
            _lines = _f.readlines()
            _headers = _lines[0].split('\t')
            _x_index, _y_index, _z_index = _headers.index('X'), _headers.index('Y'), _headers.index('Z')
            for _line in _lines[1:]:
                _line = _line.split('\t')
                _x, _y, _z = float(_line[_x_index]), float(_line[_y_index]), float(_line[_z_index])
                _cells.append((_x, _y, _z))
    finally:
        _f.close()

    return _cells


def objects_series_file_data(_experiment, _series):
    _objects_by_time = list([None] * len(paths.text_files(paths.objects(_experiment, _series))))
    for _time_point_file in paths.text_files(paths.objects(_experiment, _series)):
        _time_point = int(str(_time_point_file.split('tp_')[1]).split('.')[0]) - 1
        _objects_by_time[_time_point] = objects_time_point_file_data(_experiment, _series, _time_point_file)

    return _objects_by_time


def objects_experiment_file_data(_experiment):
    return {_series: objects_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.objects(_experiment))}


def objects_z_time_point_file_data(_experiment, _series, _group, _time_point):
    _file_path = paths.objects_z(_experiment, _series, _group, _time_point)
    try:
        _cells = []
        with open(_file_path) as _f:
            _lines = _f.readlines()
            _headers = _lines[0].split('\t')
            _x_index, _y_index, _z_index = _headers.index('X'), _headers.index('Y'), _headers.index('Z')
            for _line in _lines[1:]:
                _line = _line.split('\t')
                _x, _y, _z = float(_line[_x_index]), float(_line[_y_index]), float(_line[_z_index])
                _cells.append((_x, _y, _z))
    finally:
        _f.close()

    return _cells


def objects_z_group_file_data(_experiment, _series, _group):
    _objects_by_time = list([None] * len(paths.text_files(paths.objects_z(_experiment, _series, _group))))
    for _time_point_file in paths.text_files(paths.objects_z(_experiment, _series, _group)):
        _time_point = int(str(_time_point_file.split('tp_')[1]).split('.')[0]) - 1
        _objects_by_time[_time_point] = objects_z_time_point_file_data(_experiment, _series, _group, _time_point_file)

    return _objects_by_time


def objects_z_series_file_data(_experiment, _series):
    return {_group: objects_z_group_file_data(_experiment, _series, _group) for
            _group in paths.folders(paths.objects_z(_experiment, _series))}


def objects_z_experiment_file_data(_experiment):
    return {_series: objects_z_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.objects_z(_experiment))}


def cell_coordinates_tracked_series_file_data(_experiment, _series):
    _file_path = paths.cell_coordinates_tracked(_experiment, _series)
    try:
        _cells = []
        with open(_file_path) as _f:
            _lines = _f.readlines()
            for _line in _lines[1:]:
                _line = _line.split('\t')
                _coordinates_by_time = []
                for _coordinates in _line:
                    if _coordinates == 'None':
                        _coordinates_by_time.append(None)
                    else:
                        _coordinates_split = _coordinates.split()
                        _x, _y, _z = float(_coordinates_split[0]), float(_coordinates_split[1]), float(_coordinates_split[2])
                        _coordinates_by_time.append((_x, _y, _z))
                _cells.append(_coordinates_by_time)
    finally:
        _f.close()

    return _cells


def cell_coordinates_tracked_experiment_file_data(_experiment):
    return {_series: cell_coordinates_tracked_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.cell_coordinates_tracked(_experiment))}


def cell_coordinates_tracked_z_group_file_data(_experiment, _series, _group):
    _file_path = paths.cell_coordinates_tracked_z(_experiment, _series, _group)
    try:
        _cells = []
        with open(_file_path) as _f:
            _lines = _f.readlines()
            for _line in _lines[1:]:
                _line = _line.split('\t')
                _coordinates_by_time = []
                for _coordinates in _line:
                    if _coordinates == 'None':
                        _coordinates_by_time.append(None)
                    else:
                        _coordinates_split = _coordinates.split()
                        _x, _y, _z = float(_coordinates_split[0]), float(_coordinates_split[1]), float(_coordinates_split[2])
                        _coordinates_by_time.append((_x, _y, _z))
                _cells.append(_coordinates_by_time)
    finally:
        _f.close()

    return _cells


def cell_coordinates_tracked_z_series_file_data(_experiment, _series):
    return {_group: cell_coordinates_tracked_z_group_file_data(_experiment, _series, _group) for
            _group in paths.folders(paths.cell_coordinates_tracked_z(_experiment, _series))}


def cell_coordinates_tracked_z_experiment_file_data(_experiment):
    return {_series: cell_coordinates_tracked_z_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.cell_coordinates_tracked_z(_experiment))}


def fibers_density_time_point_file_data(_experiment, _series, _group, _z_group, _time_point):
    _file_path = paths.fibers_density(_experiment, _series, _group, _z_group, _time_point)
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


def fibers_density_z_group_file_data(_experiment, _series, _group, _z_group):
    return {_time_point: fibers_density_time_point_file_data(_experiment, _series, _group, _z_group, _time_point) for
            _time_point in paths.text_files(paths.fibers_density(_experiment, _series, _group, _z_group))}


def fibers_density_group_file_data(_experiment, _series, _group):
    return {_z_group: fibers_density_z_group_file_data(_experiment, _series, _group, _z_group) for
            _z_group in paths.folders(paths.fibers_density(_experiment, _series, _group))}


def fibers_density_series_file_data(_experiment, _series):
    return {_group: fibers_density_group_file_data(_experiment, _series, _group) for
            _group in paths.folders(paths.fibers_density(_experiment, _series))}


def fibers_density_experiment_file_data(_experiment):
    return {_series: fibers_density_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.fibers_density(_experiment))}


def normalization_file_data(_experiment, _series):
    # TODO: create it
    return False


def objects_file_names(_experiment, _series):
    _experiment_path = paths.objects(_experiment)
    _series_path = os.path.join(_experiment_path)


def serieses(_experiment_path):
    return [_series for _series in os.listdir(_experiment_path) if
            os.path.isdir(os.path.join(_experiment_path, _series)) and _series.startswith('Series')]


def group_normalization(_series_path, _group):
    _group_path = os.path.join(_series_path, _group)
    return fibers_density_time_point_file_data(_group_path)[0]


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
