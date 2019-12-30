import math

import numpy as np

from configurations import CELL_DIAMETER, EPSILON


def correlation(_array1, _array2):
    return np.corrcoef(_array1, _array2)[0][1]


def derivative(_array, _n=1):
    return list(np.diff(_array, n=_n))


def cells_distance(_properties):
    _c1_x = _properties['time_points'][0]['left_cell']['coordinates']['x']
    _c2_x = _properties['time_points'][0]['right_cell']['coordinates']['x']

    return round(abs(_c2_x - _c1_x) / CELL_DIAMETER, 10)


def is_in_region(_x, _y, _roi):
    return _roi[0] - EPSILON <= _x <= _roi[2] + EPSILON and _roi[1] - EPSILON <= _y <= _roi[3] + EPSILON


def two_lines_intersection(_l1, _l2):
    _d = _l1[0] * _l2[1] - _l1[1] * _l2[0]
    _dx = _l1[2] * _l2[1] - _l1[1] * _l2[2]
    _dy = _l1[0] * _l2[2] - _l1[2] * _l2[0]

    if _d != 0:
        x = _dx / _d
        y = _dy / _d
        return round(x, 10), round(y, 10)
    else:
        return False


def is_value_between_values(_v1, _v2, _v):
    return _v1 - EPSILON <= _v <= _v2 + EPSILON or _v2 - EPSILON <= _v <= _v1 + EPSILON


def roi(_length_x, _length_y, _offset_x, _offset_y, _cell_coordinates, _cell_diameter, _direction):
    if _direction == 'right':
        _x1 = round(_cell_coordinates['x'] + _cell_diameter / 2 + _offset_x, 10)
    elif _direction == 'left':
        _x1 = round(_cell_coordinates['x'] - _cell_diameter / 2 - _length_x, 10)
    else:
        raise Exception('No such direction ' + _direction)
    _x2 = round(_x1 + _length_x, 10)
    _y1 = round(_cell_coordinates['y'] - _length_y / 2 + _offset_y, 10)
    _y2 = round(_y1 + _length_y, 10)

    return _x1, _y1, _x2, _y2


def line(_x1, _y1, _x2, _y2):
    return round(_y1 - _y2, 10), round(_x2 - _x1, 10), - round(_x1 * _y2 - _x2 * _y1, 10)


def fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _roi, _roi_lines):
    # check if in region
    _i1_in = is_in_region(_x1, _y1, _roi)
    _i2_in = is_in_region(_x2, _y2, _roi)

    # both in region
    if _i1_in and _i2_in:
        return _original_fiber_length
    else:
        # fiber line
        _f = line(_x1, _y1, _x2, _y2)

        # get intersections
        _is = [two_lines_intersection(_f, _roi_l) for _roi_l in _roi_lines]

        # add only not false and valid
        _intersections = [
            _i for _i in _is if _i and
                                is_value_between_values(_x1, _x2, _i[0]) and
                                is_value_between_values(_y1, _y2, _i[1]) and
                                is_in_region(_i[0], _i[1], _roi)
        ]

        # remove duplicates
        _intersections = list(set(_intersections))

        # one in region
        if len(_intersections) == 1:
            _px, _py = _intersections[0]
            _inside_length = math.hypot(_px - _x1, _py - _y1)\
                if is_in_region(_x1, _y1, _roi) else math.hypot(_px - _x2, _py - _y2)
            return _inside_length * (_original_fiber_length / _current_fiber_length)

        # both out region and valid
        elif len(_intersections) == 2:
            _inside_length = math.hypot(_intersections[1][0] - _intersections[0][0],
                                        _intersections[1][1] - _intersections[0][1])
            return _inside_length * (_original_fiber_length / _current_fiber_length)

    return 0
