import math

from libs.compute_lib import two_lines_intersection, is_value_between_values, window
from libs.simulations.config import CELL_DIAMETER, EPSILON, NET_DIMENSIONS
from libs.simulations.load import elements, fiber_lengths, intersections, properties


def quantification_window_fiber_density(_simulation, _time_point, _window):
    _elements = elements(_simulation)
    _time_point_0_fiber_lengths = fiber_lengths(_simulation, _time_point=0)
    _time_point_fiber_lengths = fiber_lengths(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point)

    _master_window = get_master_window([_window])
    _elements, _time_point_0_fiber_lengths, _time_point_fiber_lengths = elements_and_intersections_by_windows(
        _elements,
        _intersections,
        _time_point_0_fiber_lengths,
        _time_point_fiber_lengths,
        _master_window
    )

    _window_lines = [
        line(_window[0], _window[1], _window[0], _window[3]),
        line(_window[0], _window[3], _window[2], _window[3]),
        line(_window[2], _window[1], _window[2], _window[3]),
        line(_window[0], _window[1], _window[2], _window[1])
    ]
    _window_size = (_window[2] - _window[0]) * (_window[3] - _window[1])
    _window_fiber_sum = 0
    for _element_index, _element in enumerate(_elements):
        _intersection_1_id, _intersection_2_id = _element
        _x1, _y1 = _intersections[_intersection_1_id]
        _x2, _y2 = _intersections[_intersection_2_id]
        _original_fiber_length = _time_point_0_fiber_lengths[_element_index]
        _current_fiber_length = _time_point_fiber_lengths[_element_index]
        _fiber = fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _window, _window_lines)
        if _fiber is not None:
            _window_fiber_sum += _fiber

    return _window_fiber_sum / _window_size


def quantification_windows_fiber_densities(_simulation, _time_point, _windows):
    _elements = elements(_simulation)
    _time_point_0_fiber_lengths = fiber_lengths(_simulation, _time_point=0)
    _time_point_fiber_lengths = fiber_lengths(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point=_time_point)

    _windows_fiber_sums = {}
    _windows_lines = {}
    _windows_sizes = {}
    _windows_to_compute = []
    for _window in _windows:
        _windows_fiber_sums[_window] = 0
        _windows_to_compute.append(_window)
        _windows_lines[_window] = [
            line(_window[0], _window[1], _window[0], _window[3]),
            line(_window[0], _window[3], _window[2], _window[3]),
            line(_window[2], _window[1], _window[2], _window[3]),
            line(_window[0], _window[1], _window[2], _window[1])
        ]
        _windows_sizes[_window] = (_window[2] - _window[0]) * (_window[3] - _window[1])

    _master_window = get_master_window(_windows)
    _elements, _time_point_0_fiber_lengths, _time_point_fiber_lengths = elements_and_intersections_by_windows(
        _elements,
        _intersections,
        _time_point_0_fiber_lengths,
        _time_point_fiber_lengths,
        _master_window
    )

    _elements_len = str(len(_elements))
    _windows_to_compute_len = str(len(_windows_to_compute))
    for _element_index, _element in enumerate(_elements):
        _intersection_1_id, _intersection_2_id = _element
        _x1, _y1 = _intersections[_intersection_1_id]
        _x2, _y2 = _intersections[_intersection_2_id]
        _original_fiber_length = _time_point_0_fiber_lengths[_element_index]
        _current_fiber_length = _time_point_fiber_lengths[_element_index]
        for _window_index, _window in enumerate(_windows_to_compute):
            _window_lines = _windows_lines[_window]
            _fiber = fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _window, _window_lines)
            if _fiber is not None:
                _windows_fiber_sums[_window] += _fiber

    for _window in _windows_to_compute:
        _windows_fiber_sums[_window] = _windows_fiber_sums[_window] / _windows_sizes[_window]

    return _windows_fiber_sums


def window_fiber_density_time_point(_arguments):
    if 'properties' not in _arguments:
        _arguments['properties'] = properties(_arguments['simulation'])
    _time_point_properties = _arguments['properties']['time_points'][_arguments['time_point']]
    if _arguments['direction'] == 'inside':
        _new_direction = 'right' if _arguments['cell_id'] == 'left_cell' else 'left'
    elif _arguments['direction'] == 'outside':
        _new_direction = 'left' if _arguments['cell_id'] == 'left_cell' else 'right'
    else:
        _new_direction = _arguments['direction']
    _time_point_window = window(
        _length_x=_arguments['length_x'] * CELL_DIAMETER,
        _length_y=_arguments['length_y'] * CELL_DIAMETER,
        _offset_x=_arguments['offset_x'] * CELL_DIAMETER,
        _offset_y=_arguments['offset_y'] * CELL_DIAMETER,
        _cell_coordinates=_time_point_properties[_arguments['cell_id']]['coordinates'],
        _cell_diameter=_time_point_properties[_arguments['cell_id']]['diameter'],
        _direction=_new_direction
    )
    if 'print' in _arguments and _arguments['print']:
        print('Computing:', _arguments['simulation'], _arguments['cell_id'], 'window', _time_point_window, 'direction',
              _arguments['direction'], 'tp', _arguments['time_point'])
    _window_fiber_density = \
        quantification_window_fiber_density(_arguments['simulation'], _arguments['time_point'], _time_point_window)

    return _arguments, _window_fiber_density


def window_fiber_density_by_time(_arguments):
    _simulation_properties = properties(_arguments['simulation'])
    _windows = []
    for _time_point in range(min(_arguments['time_points'], len(_simulation_properties['time_points']))):
        _arguments['time_point'] = _time_point
        _windows.append(window_fiber_density_time_point(_arguments)[1])

    return _arguments, _windows


def pair_distance(_properties):
    _c1_x = _properties['time_points'][0]['left_cell']['coordinates']['x']
    _c2_x = _properties['time_points'][0]['right_cell']['coordinates']['x']

    return round(abs(_c2_x - _c1_x) / CELL_DIAMETER, 10)


def is_in_window(_x, _y, _window):
    return _window[0] - EPSILON <= _x <= _window[2] + EPSILON and _window[1] - EPSILON <= _y <= _window[3] + EPSILON


def line(_x1, _y1, _x2, _y2):
    return round(_y1 - _y2, 10), round(_x2 - _x1, 10), - round(_x1 * _y2 - _x2 * _y1, 10)


def fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _window, _window_lines):
    # check if in window
    _intersection_1_in = is_in_window(_x1, _y1, _window)
    _intersection_2_in = is_in_window(_x2, _y2, _window)

    # both in window
    if _intersection_1_in and _intersection_2_in:
        return _original_fiber_length
    else:
        # fiber line
        _fiber = line(_x1, _y1, _x2, _y2)

        # get intersections
        _intersections = [two_lines_intersection(_fiber, _window_line) for _window_line in _window_lines]

        # add only not false and valid
        _intersections = [
            _intersection for _intersection in _intersections if _intersection and
            is_value_between_values(_x1, _x2, _intersection[0]) and
            is_value_between_values(_y1, _y2, _intersection[1]) and
            is_in_window(_intersection[0], _intersection[1], _window)
        ]

        # remove duplicates
        _intersections = list(set(_intersections))

        # one in window
        if len(_intersections) == 1:
            _px, _py = _intersections[0]
            _inside_length = math.hypot(_px - _x1, _py - _y1) \
                if is_in_window(_x1, _y1, _window) else math.hypot(_px - _x2, _py - _y2)
            return _inside_length * (_original_fiber_length / _current_fiber_length)

        # both out window and valid
        elif len(_intersections) == 2:
            _inside_length = math.hypot(_intersections[1][0] - _intersections[0][0],
                                        _intersections[1][1] - _intersections[0][1])
            return _inside_length * (_original_fiber_length / _current_fiber_length)

    return None


def get_expanded_window(_window):
    _expand_factor = 0.1
    return [
        _window[0] - _expand_factor,
        _window[1] - _expand_factor,
        _window[2] + _expand_factor,
        _window[3] + _expand_factor
    ]


def get_master_window(_windows):
    _x1, _y1, _x2, _y2 = -NET_DIMENSIONS['x1'], -NET_DIMENSIONS['y1'], -NET_DIMENSIONS['x2'], -NET_DIMENSIONS['y2']
    for _window in _windows:
        _x1, _y1, _x2, _y2 = min(_x1, _window[0]), min(_y1, _window[1]), max(_x2, _window[2]), max(_y2, _window[3])

    return get_expanded_window([_x1, _y1, _x2, _y2])


def elements_and_intersections_by_windows(_elements, _intersections, _time_point_0_fiber_lengths, _time_point_fiber_lengths, _window):
    _filtered_elements, _filtered_time_point_0_fiber_lengths, _filtered_time_point_fiber_lengths = [], [], []
    for _element_index, _element in enumerate(_elements):
        _intersection_1_id, _intersection_2_id = _element
        _x1, _y1 = _intersections[_intersection_1_id]
        _x2, _y2 = _intersections[_intersection_2_id]
        # both intersections in window
        if is_in_window(_x1, _y1, _window) and is_in_window(_x2, _y2, _window):
            _filtered_elements.append(_element)
            if _time_point_0_fiber_lengths is not None:
                _filtered_time_point_0_fiber_lengths.append(_time_point_0_fiber_lengths[_element_index])
            if _time_point_fiber_lengths is not None:
                _filtered_time_point_fiber_lengths.append(_time_point_fiber_lengths[_element_index])

    return _filtered_elements, _filtered_time_point_0_fiber_lengths, _filtered_time_point_fiber_lengths
