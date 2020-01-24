import gzip
import json
import os

import pickle

from libs.experiments import paths


def time_point_file_name_to_number(_time_point_file_name):
    return int(str(_time_point_file_name.split('tp_')[1]).split('.')[0])


def series_file_name_to_name(_series_file_name):
    return 'Series ' + str(str(_series_file_name.split('series_')[1]).split('.')[0])


def experiment_groups_as_tuples(_experiment):
    _tuples = []
    for _series in paths.folders(paths.structured(_experiment)):
        for _group in paths.folders(paths.structured(_experiment, _series)):
            _tuples.append((_experiment, int(_series.split()[1]), _group))

    return _tuples


def experiments_groups_as_tuples(_experiments):
    _tuples = []
    for _experiment in _experiments:
        _tuples += experiment_groups_as_tuples(_experiment)

    return _tuples


def image_properties(_experiment, _series):
    _file_path = paths.image_properties(_experiment, 'series_' + str(_series.split()[1]) + '.json')
    try:
        with open(_file_path) as _json:
            return json.load(_json)
    finally:
        _json.close()


def group_properties(_experiment, _series_id, _group):
    _file_path = paths.group_properties(_experiment, 'Series ' + str(_series_id), _group)
    try:
        with open(_file_path) as _json:
            return json.load(_json)
    finally:
        _json.close()


def structured_image(_experiment, _series_id, _group, _time_point):
    _image_path = paths.structured(_experiment, 'Series ' + str(_series_id), _group, str(_time_point) + '.pkl')
    try:
        with gzip.open(_image_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


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
        _time_point = time_point_file_name_to_number(_time_point_file)
        _objects_by_time[_time_point - 1] = objects_time_point_file_data(_experiment, _series, _time_point_file)

    return _objects_by_time


def objects_experiment_file_data(_experiment):
    return {_series: objects_series_file_data(_experiment, _series) for
            _series in paths.folders(paths.objects(_experiment))}


def cell_coordinates_tracked_series_file_data(_experiment, _series):
    _file_path = paths.cell_coordinates_tracked(_experiment, _series)
    try:
        _cells = []
        with open(_file_path) as _f:
            _lines = _f.readlines()
            for _line in _lines:
                _line = _line.replace('\n', '')
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


def fibers_densities(_experiment, _series_id, _group, _time_point):
    _fibers_densities_path = paths.fibers_densities(_experiment, 'Series ' + str(_series_id), _group, str(_time_point) + '.pkl')
    if os.path.isfile(_fibers_densities_path):
        if os.path.getsize(_fibers_densities_path) == 0:
            return {}
        try:
            with gzip.open(_fibers_densities_path, 'rb') as _pickle:
                return pickle.load(_pickle)
        finally:
            _pickle.close()
    else:
        return {}


def normalization_series_file_data(_experiment, _series):
    _file_path = paths.normalization(_experiment, _series)
    try:
        with open(_file_path, 'r') as _file:
            _lines = _file.readlines()
            _line = _lines[0].split()
            _average, _std = float(_line[0]), float(_line[1])
    finally:
        _file.close()

    return _average, _std


def normalization_experiment_file_data(_experiment):
    return {series_file_name_to_name(_series):
            normalization_series_file_data(_experiment, series_file_name_to_name(_series)) for
            _series in paths.text_files(paths.normalization(_experiment))}
