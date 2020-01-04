import numpy as np

from libs.simulations.config import EPSILON


def correlation(_array1, _array2):
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
