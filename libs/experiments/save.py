import json
import os

from libs import save_lib
from libs.experiments import paths
from libs.save_lib import to_pickle


def cell_coordinates_tracked(_experiment, _series_id, _cell_coordinates):
    _experiment_path = paths.cell_coordinates_tracked(_experiment)
    os.mkdir(_experiment_path) if not os.path.isdir(_experiment_path) else None
    _file_path = os.path.join(_experiment_path, 'series_' + str(_series_id) + '.txt')
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


def normalization_line(_experiment, _series, _line):
    _experiment_path = paths.normalization_lines(_experiment)
    os.mkdir(_experiment_path) if not os.path.isdir(_experiment_path) else None
    _file_path = os.path.join(_experiment_path, 'series_' + str(_series.split()[1]) + '.txt')
    try:
        with open(_file_path, 'w') as _f:
            _f.write(str(_line[0]) + '\t' + str(_line[1]) + '\t' + str(_line[2]) + '\t' + str(_line[3]))
    finally:
        _f.close()


def image_properties(_experiment, _series_id, _image_properties):
    _experiment_path = paths.image_properties(_experiment)
    os.mkdir(_experiment_path) if not os.path.isdir(_experiment_path) else None
    _path = os.path.join(_experiment_path, 'series_' + str(_series_id) + '.json')
    save_lib.to_json(_image_properties, _path)


def fibers_densities(_experiment, _series_id, _group, _time_point, _fibers_densities):
    _fibers_densities_path = paths.fibers_densities(_experiment, 'Series ' + str(_series_id), _group, str(_time_point) + '.pkl')
    to_pickle(_fibers_densities, _fibers_densities_path)


def blacklist(_experiment, _series_id, _group, _blacklist):
    _path = paths.blacklist(_experiment, _series_id, _group)
    to_pickle(_blacklist, _path)
