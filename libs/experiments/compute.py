import math

import numpy as np
from scipy.ndimage import rotate

from libs.compute_lib import roi, z_score
from libs.experiments import load, save
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def z_score_fibers_density_array(_fibers_density, _normalization):
    _average, _std = _normalization
    return {k: [z_score(_x, _average, _std) for _x in _fibers_density[k]] for k in _fibers_density.keys()}


def fibers_density_cut_edges(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:-_cut_amount]) for k in _fibers_density.keys()}


def fibers_density_cut_left_edge(_fibers_density, _cut_amount=4):
    return {k: list(_fibers_density[k][_cut_amount:]) for k in _fibers_density.keys()}


def cells_distance_in_cell_size(_experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates):
    _image_properties = load.image_properties(_experiment, 'Series ' + str(_series_id))
    _image_resolutions = _image_properties['resolutions']
    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
    _x1, _y1, _z1 = _x1 * _image_resolutions['x'], _y1 * _image_resolutions['y'], _z1 * _image_resolutions['z']
    _x2, _y2, _z2 = _x2 * _image_resolutions['x'], _y2 * _image_resolutions['y'], _z2 * _image_resolutions['z']

    return math.sqrt((_x1 - _x2) ** 2 + (_y1 - _y2) ** 2 + (_z1 - _z2) ** 2) / CELL_DIAMETER_IN_MICRONS


def angle_between_three_points(_a, _b, _c):
    _angle = abs(math.degrees(math.atan2(_c[1] - _b[1], _c[0] - _b[0]) - math.atan2(_a[1] - _b[1], _a[0] - _b[0])))
    return _angle if _b[1] <= _a[1] else -_angle


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


def axes_padding(_2d_image_shape, _angle):
    _image_zeros = np.zeros(_2d_image_shape)
    _image_zeros_rotated = rotate(_image_zeros, _angle)
    _image_zeros_shape = _image_zeros_rotated.shape

    # return x, y
    return int(round((_image_zeros_shape[1] - _2d_image_shape[1]) / 2)), \
        int(round((_image_zeros_shape[0] - _2d_image_shape[0]) / 2))


def image_center_coordinates(_image_shape):
    return [_value / 2 for _value in _image_shape]


def roi_by_microns(_resolution_x, _resolution_y, _resolution_z, _length_x, _length_y, _length_z, _offset_x, _offset_y,
                   _offset_z, _cell_coordinates, _direction):
    _cell_diameter_in_pixels = CELL_DIAMETER_IN_MICRONS / _resolution_x
    _length_x_in_pixels = _length_x / _resolution_x
    _length_y_in_pixels = _length_y / _resolution_y
    _length_z_in_pixels = _length_z / _resolution_z
    _offset_x_in_pixels = _offset_x / _resolution_x
    _offset_y_in_pixels = _offset_y / _resolution_y
    _offset_z_in_pixels = _offset_z / _resolution_z

    _x1, _y1, _x2, _y2 = [int(round(_value)) for _value in roi(
        _length_x=_length_x_in_pixels,
        _length_y=_length_y_in_pixels,
        _offset_x=_offset_x_in_pixels,
        _offset_y=_offset_y_in_pixels,
        _cell_coordinates={
            'x': _cell_coordinates['x'],
            'y': _cell_coordinates['y']
        },
        _cell_diameter=_cell_diameter_in_pixels,
        _direction=_direction
    )]

    _z1 = int(round(_cell_coordinates['z'] - _length_z_in_pixels / 2 + _offset_z, 10))
    _z2 = int(round(_z1 + _length_z_in_pixels, 10))

    return _x1, _y1, _z1, _x2, _y2, _z2


def roi_fibers_density(_experiment, _series, _group, _time_point, _roi):
    _time_point_image = load.structured_image(_experiment, _series, _group, _time_point)
    _x1, _y1, _z1, _x2, _y2, _z2 = _roi
    return np.mean(_time_point_image[_z1:_z2, _y1:_y2, _x1:_x2])


def roi_fibers_density_time_point(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x,
                                  _offset_y, _offset_z, _cell_id, _direction, _time_point, _group_properties=None):
    _group_properties = _group_properties if _group_properties is not None else load.group_properties(_experiment,
                                                                                                      _series_id,
                                                                                                      _group)
    _time_point_properties = _group_properties['time_points'][_time_point]
    _time_point_fibers_densities = load.fibers_densities(_experiment, _series_id, _group, _time_point)
    _time_point_roi = roi_by_microns(
        _resolution_x=_group_properties['time_points'][_time_point]['resolutions']['x'],
        _resolution_y=_group_properties['time_points'][_time_point]['resolutions']['y'],
        _resolution_z=_group_properties['time_points'][_time_point]['resolutions']['z'],
        _length_x=_length_x,
        _length_y=_length_y,
        _length_z=_length_z,
        _offset_x=_offset_x,
        _offset_y=_offset_y,
        _offset_z=_offset_z,
        _cell_coordinates=_group_properties['time_points'][_time_point][_cell_id]['coordinates'],
        _direction='right' if
        (_cell_id, _direction) == ('left_cell', 'inside') or
        (_cell_id, _direction) == ('right_cell', 'outside') or
        (_cell_id, _direction) == ('cell', 'right') else 'left'
    )
    if _time_point_roi in _time_point_fibers_densities:
        return _time_point_fibers_densities[_time_point_roi]
    else:
        print('Computing:', _experiment, _series_id, _group, _cell_id, 'roi', _time_point_roi, 'direction', _direction,
              'tp', _time_point)
        _roi_fibers_density = roi_fibers_density(_experiment, _series_id, _group, _time_point, _time_point_roi)
        _time_point_fibers_densities[_time_point_roi] = _roi_fibers_density
        save.fibers_densities(_experiment, _series_id, _group, _time_point, _time_point_fibers_densities)
        return _roi_fibers_density


def roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x, _offset_y,
                               _offset_z, _cell_id, _direction, _time_points):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    return [
        roi_fibers_density_time_point(
            _experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x, _offset_y, _offset_z, _cell_id,
            _direction, _time_point, _group_properties
        ) for _time_point in range(min(_time_points, len(_group_properties['time_points'])))
    ]


def roi_fibers_density_by_time_pairs(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x,
                                     _offset_y, _offset_z, _direction, _time_points):
    return {
        'left_cell': roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z,
                                                _offset_x, _offset_y, _offset_z, 'left_cell', _direction, _time_points),
        'right_cell': roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z,
                                                 _offset_x, _offset_y, _offset_z, 'right_cell', _direction,
                                                 _time_points),
    }


def minimum_time_points(_experiments_tuples):
    _minimum_time_points = 10000
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _minimum_time_points = min(_minimum_time_points, len(_group_properties['time_points']))

    return _minimum_time_points
