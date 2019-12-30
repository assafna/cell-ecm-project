import math
import os

import numpy as np

from libs.experiments import config, paths

MAX_DISTANCE_CHANGE = 20


def load_objects(_path):
    _files = [_file for _file in os.listdir(_path) if _file.startswith('tp_') and _file.endswith('.txt')]
    _files_data = list([None] * len(_files))
    for _file_name in _files:
        _file_path = os.path.join(_path, _file_name)
        _tp = int(str(_file_name.split('tp_')[1]).split('.')[0]) - 1
        _cells = []
        try:
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
        _files_data[_tp] = _cells

    return _files_data


def load_serieses(_experiment):
    _path = os.path.join(paths.MANIPULATIONS, _experiment)
    return [_series for _series in os.listdir(_path) if
            os.path.isdir(os.path.join(_path, _series)) and _series.startswith('Series')]


def load_z_groups(_experiment, _series):
    _path = os.path.join(paths.objects(_experiment, True), _series)
    return [_group for _group in os.listdir(_path) if
            os.path.isdir(os.path.join(_path, _group)) and _group.startswith('cell')]


def load_objects_data_by_time(_experiment, _is_z):
    _objects_data = {}
    for _series in load_serieses(_experiment):
        _series_path = os.path.join(paths.objects(_experiment, _is_z), _series)
        if not os.path.isdir(_series_path):
            print('Series not found!', _series)
            return None
        print(_series)
        if _is_z:
            _z_group_objects_data = {}
            for _z_group in load_z_groups(_experiment, _series):
                print(_z_group)
                _z_group_path = os.path.join(_series_path, _z_group)
                _z_group_objects_data[_z_group] = load_objects(_z_group_path)
            _objects_data[_series] = _z_group_objects_data
        else:
            _objects_data[_series] = load_objects(_series_path)

    return _objects_data


def cell_coordinates(_objects_data):
    _cell_coordinates = [[_coordinates] for _coordinates in _objects_data[0]]
    for _tp in _objects_data[1:]:
        for _cell_index, _cell in enumerate(_cell_coordinates):
            if _cell[-1] is not None:
                _previous_cell_x, _previous_cell_y, _previous_cell_z = _cell[-1]
                _distances = [math.sqrt(
                    math.pow(_previous_cell_x - _optional_cell_x, 2) +
                    math.pow(_previous_cell_y - _optional_cell_y, 2) +
                    math.pow(_previous_cell_z - _optional_cell_z, 2)
                ) for _optional_cell_x, _optional_cell_y, _optional_cell_z in _tp]
                _min_distance_index = int(np.argmin(_distances)) if len(_distances) > 0 else None
                _cell_coordinates[_cell_index].append(
                    _tp[_min_distance_index] if _min_distance_index is not None and min(
                        _distances) < MAX_DISTANCE_CHANGE else None
                )
            else:
                _cell_coordinates[_cell_index].append(None)

    return _cell_coordinates


def cell_coordinates_by_time(_objects_data, _is_z):
    _cell_coordinates = {}
    for _series in _objects_data:
        _series_cell_coordinates = {}
        _series_objects_data = _objects_data[_series]
        if _is_z:
            for _z_group in _series_objects_data:
                _series_cell_coordinates[_z_group] = cell_coordinates(_series_objects_data[_z_group])
        else:
            _series_cell_coordinates = cell_coordinates(_series_objects_data)
        _cell_coordinates[_series] = _series_cell_coordinates

    return _cell_coordinates


def save(_cell_coordinates, _file_path):
    _lines = ''
    for _cell_index, _cell in enumerate(_cell_coordinates):
        _line = ''
        for _tp_index, _tp in enumerate(_cell):
            _tp_str = str(_tp[0]) + ' ' + str(_tp[1]) + ' ' + str(_tp[2]) if _tp is not None else 'None'
            _line += _tp_str + '\t' if _tp_index < len(_cell) - 1 else _tp_str
        _lines += _line + '\n' if _cell_index < len(_cell_coordinates) - 1 else _line
    try:
        with open(_file_path, 'w') as _f:
            _f.write(_lines)
    finally:
        _f.close()


def save_to_file(_cell_coordinates, _experiment, _is_z):
    _experiment_path = paths.cell_coordinates_tracked(_experiment, _is_z)
    for _series in _cell_coordinates:
        _series_cell_coordinates = _cell_coordinates[_series]
        if _is_z:
            for _z_group in _series_cell_coordinates:
                _z_group_cell_coordinates = _series_cell_coordinates[_z_group]
                _series_path = os.path.join(_experiment_path, _series)
                os.mkdir(_series_path) if not os.path.isdir(_series_path) else None
                _file_path = os.path.join(_series_path, _z_group + '.txt')
                save(_z_group_cell_coordinates, _file_path)
        else:
            _file_path = os.path.join(_experiment_path, 'series_' + str(_series.split()[1]) + '.txt')
            save(_series_cell_coordinates, _file_path)


if __name__ == '__main__':
    _is_z = False
    for experiment in config.SINGLE_CELL:
        print('Experiment', experiment)
        os.mkdir(paths.cell_coordinates_tracked(experiment, _is_z)) if not os.path.isdir(
            paths.cell_coordinates_tracked(experiment, _is_z)) else None
        objects = load_objects_data_by_time(experiment, _is_z=_is_z)
        coordinates = cell_coordinates_by_time(objects, _is_z=_is_z) if objects is not None else None
        save_to_file(coordinates, experiment, _is_z=_is_z) if coordinates is not None else None
