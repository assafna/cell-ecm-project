import math

from libs.experiments import load
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def z_score(_x, _average, _std):
    return (_x - _average) / _std


def z_score_fibers_density_array(_fibers_density, _normalization):
    _average, _std = _normalization
    return {k: [z_score(_x, _average, _std) for _x in _fibers_density[k]] for k in _fibers_density.keys()}


def fibers_density_cut_edges(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:-_cut_amount]) for k in _fibers_density.keys()}


def fibers_density_cut_left_edge(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:]) for k in _fibers_density.keys()}


def cells_distance_in_cell_size(_experiment, _series, _cell_1_coordinates, _cell_2_coordinates):
    _series_id = _series.split()[1]
    _image_properties = load.image_properties(_experiment, _series)
    _image_resolutions = _image_properties['resolutions']
    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
    _x1, _y1, _z1 = _x1 * _image_resolutions['x'], _y1 * _image_resolutions['y'], _z1 * _image_resolutions['z']
    _x2, _y2, _z2 = _x2 * _image_resolutions['x'], _y2 * _image_resolutions['y'], _z2 * _image_resolutions['z']

    return math.sqrt((_x1 - _x2) ** 2 + (_y1 - _y2) ** 2 + (_z1 - _z2) ** 2) / CELL_DIAMETER_IN_MICRONS


def angle_between_three_points(_a, _b, _c):
    return round(abs(math.degrees(math.atan2(_c[1] - _b[1], _c[0] - _b[0]) - math.atan2(_a[1] - _b[1], _a[0] - _b[0]))))


def rotate_point_around_another_point(_point, _angle_in_radians, _around_point):
    _z = None
    if len(_point) == 3:
        _x, _y, _z = _point
    else:
        _x, _y = _point
    _offset_x, _offset_y = _around_point
    _adjusted_x = (_x - _offset_x)
    _adjusted_y = (_y - _offset_y)
    _cos_rad = math.cos(_angle_in_radians)
    _sin_rad = math.sin(_angle_in_radians)
    _qx = int(round(_offset_x + _cos_rad * _adjusted_x + _sin_rad * _adjusted_y))
    _qy = int(round(_offset_y + -_sin_rad * _adjusted_x + _cos_rad * _adjusted_y))

    return [_qx, _qy, _z] if len(_point) == 3 else [_qx, _qy]
