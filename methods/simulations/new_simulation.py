import math
import os
import time

import libs.simulations.load
from libs.simulations import load, paths
from libs.simulations.config import CELL_DIAMETER, NET_DIMENSIONS, ORIGIN_COORDINATES, \
    ELEMENTS_FILE_NAME
from libs.simulations.save import to_pickle


def create_properties(_simulation):
    print(_simulation, 'Creating properties')
    _simulation_path = paths.raw(_simulation)
    _cells_diameters = load.cells_diameters(_simulation)
    _time_points = load.time_points(_simulation)
    _time_points_data = []
    if 'D_' in _simulation:
        _distance = float(_simulation.split('D')[0])
        _left_cell_x = round(-CELL_DIAMETER / 2 * _distance, 10)
        _right_cell_x = round(CELL_DIAMETER / 2 * _distance, 10)

        for _time_point in _time_points:
            _left_cell_diameter = _cells_diameters[_time_point][0]
            _right_cell_diameter = _cells_diameters[_time_point][1]
            _time_point_line = {
                "left_cell": {"coordinates": {"x": _left_cell_x, "y": 0.0}, "diameter": _left_cell_diameter},
                "right_cell": {"coordinates": {"x": _right_cell_x, "y": 0.0}, "diameter": _right_cell_diameter}
            }
            _time_points_data.append(_time_point_line)
    else:
        for _time_point in _time_points:
            _cell_diameter = _cells_diameters[_time_point]
            _time_point_line = {
                "cell": {"coordinates": {"x": 0.0, "y": 0.0}, "diameter": _cell_diameter}
            }
            _time_points_data.append(_time_point_line)

    # prepare properties structure
    _data = {
        'name': _simulation,
        'net_dimensions': NET_DIMENSIONS,
        'origin_coordinates': ORIGIN_COORDINATES,
        'time_points': _time_points_data
    }

    return _data


def create_elements(_simulation):
    print(_simulation, 'Creating elements')
    _simulation_path = paths.raw(_simulation)
    _elements_path = os.path.join(_simulation_path, ELEMENTS_FILE_NAME)
    _elements = []
    try:
        with open(_elements_path, 'r') as _f:
            _l = _f.readline()
            while _l:
                _l_split = [_x.strip() for _x in _l.split(',')]
                _element = [int(_l_split[1]), int(_l_split[2])]
                _elements.append(_element)
                _l = _f.readline()
    finally:
        _f.close()

    return _elements


def create_time_point(_simulation, _time_point):
    _simulation_path = paths.raw(_simulation)
    _time_point_path = os.path.join(_simulation_path, str(_time_point) + '.txt')
    _intersections = [0]
    try:
        with open(_time_point_path, 'r') as _f:
            _l = _f.readline()
            while _l:
                _l_split = [_x.strip() for _x in _l.split(',')]
                _i_id, _x, _y = _l_split
                _intersections.append([float(_x), float(_y)])
                _l = _f.readline()
    finally:
        _f.close()

    return _intersections


def create_fiber_lengths(_elements, _intersections):
    _fiber_lengths = []
    for _element in _elements:
        _i1_id, _i2_id = _element
        _x1, _y1 = _intersections[_i1_id]
        _x2, _y2 = _intersections[_i2_id]
        _fiber_lengths.append(math.hypot(_x2 - _x1, _y2 - _y1))

    return _fiber_lengths


def process_simulation(_simulation, _overwrite=False):
    _start_time = time.time()

    _simulation_structured_path = paths.structured(_simulation)
    os.makedirs(_simulation_structured_path) if not os.path.isdir(_simulation_structured_path) else None

    # properties
    _properties_pickle_path = os.path.join(_simulation_structured_path, 'properties.pkl')
    if not _overwrite and not os.path.isfile(_properties_pickle_path):
        _properties = create_properties(_simulation)
        print(_simulation, 'Pickling properties')
        to_pickle(_properties, _properties_pickle_path)

    # elements
    _elements_pickle_path = os.path.join(_simulation_structured_path, 'elements.pkl')
    if not _overwrite and not os.path.isfile(_elements_pickle_path):
        _elements = create_elements(_simulation)
        print(_simulation, 'Pickling elements')
        to_pickle(_elements, _elements_pickle_path)

    # time points
    for _time_point in load.time_points(_simulation):
        _time_point_pickle_path = os.path.join(_simulation_structured_path, str(_time_point) + '.pkl')
        if not _overwrite and not os.path.isfile(_time_point_pickle_path):
            _time_point_data = create_time_point(_simulation, _time_point)
            print(_simulation, 'Pickling time-point:', _time_point)
            to_pickle(_time_point_data, _time_point_pickle_path)

    # fiber lengths
    _simulation_fiber_lengths_path = paths.fibers_lengths(_simulation)
    os.makedirs(_simulation_fiber_lengths_path) if not os.path.isdir(_simulation_fiber_lengths_path) else None
    _elements = None
    for _time_point in load.time_points(_simulation):
        _fiber_lengths_pickle_path = os.path.join(_simulation_fiber_lengths_path, str(_time_point) + '.pkl')
        if not _overwrite and not os.path.isfile(_fiber_lengths_pickle_path):
            _elements = load.elements(_simulation) if _elements is None else _elements
            _intersections = load.intersections(_simulation, _time_point)
            _time_point_fiber_lengths = create_fiber_lengths(_elements, _intersections)
            print(_simulation, 'Pickling fiber lengths time-point:', _time_point)
            to_pickle(_time_point_fiber_lengths, _fiber_lengths_pickle_path)

    print(_simulation, 'Finished!', 'Total ' + str(round(time.time() - _start_time, 2)) + ' seconds')


def process_simulations(_simulations, _overwrite=False):
    for _simulation in _simulations:
        process_simulation(_simulation, _overwrite)


def process_all_simulations(_overwrite=False):
    process_simulations(libs.simulations.load.raw(), _overwrite)


process_all_simulations()
