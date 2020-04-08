import math

import numpy as np

from libs.simulations.config import EPSILON


def correlation(_array1, _array2):
    # TODO: change to pearson correlation with p value
    return np.corrcoef(_array1, _array2)[0][1]


def derivative(_array, _n=1):
    return list(np.diff(_array, n=_n))


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
        _y1 = round(_cell_coordinates['y'] - _length_y / 2 + _offset_y, 10)
    elif _direction == 'left':
        _x1 = round(_cell_coordinates['x'] - _cell_diameter / 2 - _offset_x - _length_x, 10)
        _y1 = round(_cell_coordinates['y'] - _length_y / 2 + _offset_y, 10)
    elif _direction == 'up':
        _x1 = round(_cell_coordinates['x'] - _length_x / 2 + _offset_x, 10)
        _y1 = round(_cell_coordinates['y'] + _cell_diameter / 2 + _offset_y, 10)
    elif _direction == 'down':
        _x1 = round(_cell_coordinates['x'] - _length_x / 2 + _offset_x, 10)
        _y1 = round(_cell_coordinates['y'] - _cell_diameter / 2 - _offset_y - _length_y, 10)
    else:
        raise Exception('No such direction ' + _direction)
    _x2 = round(_x1 + _length_x, 10)
    _y2 = round(_y1 + _length_y, 10)

    return _x1, _y1, _x2, _y2


def z_score(_x, _average, _std):
    return (_x - _average) / _std


def z_score_fibers_densities_array(_fibers_densities, _normalization):
    _average, _std = _normalization
    return [z_score(_x, _average, _std) for _x in _fibers_densities]


def p_value_text(_p_value):
    if _p_value is None:
        return None

    _values = [0.0001, 0.001, 0.01, 0.05]
    for _value in _values:
        if _p_value < _value:
            return 'p-value < ' + str(_value)

    return 'p-value > ' + str(_values[-1])


def distance_from_a_point_to_a_line(_line, _point):
    _x1, _y1, _x2, _y2 = _line
    _x0, _y0 = _point

    return abs((_y2 - _y1) * _x0 - (_x2 - _x1) * _y0 + _x2 * _y1 - _y2 * _x1) / \
        math.sqrt((_y2 - _y1) ** 2 + (_x2 - _x1) ** 2)
