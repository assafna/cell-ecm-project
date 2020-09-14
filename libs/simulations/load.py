import os

from libs import load_lib
from libs.simulations import paths
from libs.simulations.config import CELLS_DIAMETERS_FILE_NAME


def raw_files(_simulation):
    _simulation_path = paths.raw(_simulation)
    return [_file for _file in os.listdir(_simulation_path) if _file.endswith('.txt')]


def pair_diameters(_simulation):
    _simulation_path = paths.raw(_simulation)
    _cells_diameters_path = os.path.join(_simulation_path, CELLS_DIAMETERS_FILE_NAME)
    _diameters = []
    with open(_cells_diameters_path, 'r') as _f:
        _l = _f.readline()
        while _l:
            if '\t' in _l:
                _l_split = [_x.strip() for _x in _l.split('\t')]
            else:
                _l_split = [_x.strip() for _x in _l.split('   ')]
            _cells_diameters = [float(_l_split[1]), float(_l_split[2])] if len(_l_split) > 2 else float(_l_split[1])
            _diameters.append(_cells_diameters)
            _l = _f.readline()

    return _diameters


def time_points(_simulation):
    _raw_files = raw_files(_simulation)
    _files = [_file for _file in _raw_files if _file.startswith('tp_') and _file.endswith('.txt')]

    return sorted([int(os.path.splitext(_file)[0].split('tp_')[1]) for _file in _files])


def properties(_simulation):
    return load_lib.from_pickle(os.path.join(paths.structured(_simulation), 'properties.pkl'))


def elements(_simulation):
    return load_lib.from_pickle(os.path.join(paths.structured(_simulation), 'elements.pkl'))


def intersections(_simulation, _time_point):
    return load_lib.from_pickle(os.path.join(paths.structured(_simulation), str(_time_point) + '.pkl'))


def fiber_lengths(_simulation, _time_point):
    return load_lib.from_pickle(os.path.join(paths.fiber_lengths(_simulation), str(_time_point) + '.pkl'))


def raw():
    return [_simulation for _simulation in os.listdir(paths.RAW)
            if os.path.isdir(os.path.join(paths.RAW, _simulation))]


def structured():
    return [_simulation for _simulation in os.listdir(paths.STRUCTURED)
            if os.path.isdir(os.path.join(paths.STRUCTURED, _simulation))]


def normalization(_simulation):
    return load_lib.from_pickle(paths.normalization(_simulation + '.pkl'))
