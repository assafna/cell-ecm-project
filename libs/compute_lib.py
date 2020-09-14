import math

import numpy as np
from scipy.stats import pearsonr

from libs.simulations.config import EPSILON


def correlation(_array_1, _array_2, _with_p_value=False):
    return pearsonr(_array_1, _array_2) if _with_p_value else pearsonr(_array_1, _array_2)[0]


def derivative(_array, _n=1):
    return list(np.diff(_array, n=_n))


def two_lines_intersection(_line_1, _line_2):
    _d = _line_1[0] * _line_2[1] - _line_1[1] * _line_2[0]
    _dx = _line_1[2] * _line_2[1] - _line_1[1] * _line_2[2]
    _dy = _line_1[0] * _line_2[2] - _line_1[2] * _line_2[0]

    if _d != 0:
        x = _dx / _d
        y = _dy / _d
        return round(x, 10), round(y, 10)
    else:
        return False


def is_value_between_values(_value_1, _value_2, _value):
    return _value_1 - EPSILON <= _value <= _value_2 + EPSILON or _value_2 - EPSILON <= _value <= _value_1 + EPSILON


def window(_length_x, _length_y, _offset_x, _offset_y, _cell_coordinates, _cell_diameter, _direction):
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


def z_score_array(_array, _average, _std):
    return [z_score(_x, _average, _std) for _x in _array]


def distance_from_a_point_to_a_line(_line, _point):
    _x1, _y1, _x2, _y2 = _line
    _x0, _y0 = _point

    return abs((_y2 - _y1) * _x0 - (_x2 - _x1) * _y0 + _x2 * _y1 - _y2 * _x1) / \
        math.sqrt((_y2 - _y1) ** 2 + (_x2 - _x1) ** 2)
