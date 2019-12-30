import gzip
import os
import pickle

from libs.simulations import paths
from libs.simulations.config import CELLS_DIAMETERS_FILE_NAME


def raw_files(_simulation):
    _simulation_path = paths.raw(_simulation)
    return [_file for _file in os.listdir(_simulation_path) if _file.endswith('.txt')]


def cells_diameters(_simulation):
    _simulation_path = paths.raw(_simulation)
    _cells_diameters_path = os.path.join(_simulation_path, CELLS_DIAMETERS_FILE_NAME)
    _diameters = []
    try:
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
    finally:
        _f.close()

    return _diameters


def time_points(_simulation):
    _raw_files = raw_files(_simulation)
    _files = [_file for _file in _raw_files if _file.startswith('tp_') and _file.endswith('.txt')]

    return sorted([int(_file.split('.')[0]) for _file in _files])


def properties(_simulation):
    _simulation_structured_path = paths.structured(_simulation)
    _properties_path = os.path.join(_simulation_structured_path, 'properties.pkl')
    try:
        with gzip.open(_properties_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


def elements(_simulation):
    _simulation_structured_path = paths.structured(_simulation)
    _elements_path = os.path.join(_simulation_structured_path, 'elements.pkl')
    try:
        with gzip.open(_elements_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


def intersections(_simulation, _time_point):
    _simulation_structured_path = paths.structured(_simulation)
    _intersections_path = os.path.join(_simulation_structured_path, str(_time_point) + '.pkl')
    try:
        with gzip.open(_intersections_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


def fibers_lengths(_simulation, _time_point):
    _simulation_fibers_lengths_path = paths.fibers_lengths(_simulation)
    _fibers_lengths_path = os.path.join(_simulation_fibers_lengths_path, str(_time_point) + '.pkl')
    try:
        with gzip.open(_fibers_lengths_path, 'rb') as _pickle:
            return pickle.load(_pickle)
    finally:
        _pickle.close()


def fibers_densities(_simulation, _time_point):
    _simulation_fibers_densities_path = paths.fibers_densities(_simulation)
    _fibers_densities_path = os.path.join(_simulation_fibers_densities_path, str(_time_point) + '.pkl')
    if os.path.isfile(_fibers_densities_path):
        try:
            with gzip.open(_fibers_densities_path, 'rb') as _pickle:
                return pickle.load(_pickle)
        finally:
            _pickle.close()
    else:
        return {}


def raw():
    return [_simulation for _simulation in os.listdir(paths.RAW)
            if os.path.isdir(os.path.join(paths.RAW, _simulation))]
