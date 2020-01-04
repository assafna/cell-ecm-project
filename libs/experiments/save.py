import os

from libs.experiments import paths


def cell_coordinates_tracked(_experiment, _series, _cell_coordinates):
    _experiment_path = paths.cell_coordinates_tracked(_experiment)
    os.mkdir(_experiment_path) if not os.path.isdir(_experiment_path) else None
    _file_path = os.path.join(_experiment_path, 'series_' + str(_series.split()[1]) + '.txt')
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


def cell_coordinates_tracked_z(_experiment, _series, _group, _cell_coordinates):
    _series_path = paths.cell_coordinates_tracked_z(_experiment, _series)
    os.mkdir(_series_path) if not os.path.isdir(_series_path) else None
    _file_path = os.path.join(_series_path, _group + '.txt')
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
