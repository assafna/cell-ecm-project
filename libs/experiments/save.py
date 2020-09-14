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
        for _time_frame_index, _time_frame in enumerate(_cell):
            _time_frame_str = str(_time_frame[0]) + ' ' + str(_time_frame[1]) + ' ' + str(_time_frame[2]) \
                if _time_frame is not None else 'None'
            _line += _time_frame_str + '\t' if _time_frame_index < len(_cell) - 1 else _time_frame_str
        _lines += _line + '\n' if _cell_index < len(_cell_coordinates) - 1 else _line
    with open(_file_path, 'w') as _f:
        _f.write(_lines)


def image_properties(_experiment, _series_id, _image_properties):
    _experiment_path = paths.image_properties(_experiment)
    os.mkdir(_experiment_path) if not os.path.isdir(_experiment_path) else None
    _path = os.path.join(_experiment_path, 'series_' + str(_series_id) + '.json')
    save_lib.to_json(_image_properties, _path)


def blacklist(_experiment, _series_id=None, _blacklist=None):
    if _blacklist is None:
        _blacklist = {}
    _path = paths.blacklist(_experiment, _series_id)
    to_pickle(_blacklist, _path)
