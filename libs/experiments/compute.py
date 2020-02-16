import math

import numpy as np
from scipy.ndimage import rotate

from libs.compute_lib import roi
from libs.experiments import load, save
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def cells_distance_in_cell_size(_experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates):
    _image_properties = load.image_properties(_experiment, _series_id)
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

    _z1 = int(round(_cell_coordinates['z'] - _length_z_in_pixels / 2 + _offset_z_in_pixels, 10))
    _z2 = int(round(_z1 + _length_z_in_pixels, 10))

    return _x1, _y1, _z1, _x2, _y2, _z2


def roi_fibers_density(_experiment, _series, _group, _time_point, _roi):
    _time_point_image = load.structured_image(_experiment, _series, _group, _time_point)
    _x1, _y1, _z1, _x2, _y2, _z2 = _roi

    _out_of_boundaries = False

    _z_shape, _y_shape, _x_shape = _time_point_image.shape
    if any([_x1 < 0, _y1 < 0, _z1 < 0, _x2 >= _x_shape, _y2 >= _y_shape, _z2 >= _z_shape]):
        _out_of_boundaries = True

        # fix boundaries
        _x1, _y1, _z1 = max(0, _x1), max(0, _y1), max(0, _z1)
        _x2, _y2, _z2 = min(_x2, _x_shape), min(_y2, _y_shape), min(_z2, _z_shape)

    _roi_pixels = _time_point_image[_z1:_z2, _y1:_y2, _x1:_x2]
    _non_zero_mask = np.nonzero(_roi_pixels)

    # check if more than 1% is black
    if not _out_of_boundaries and np.count_nonzero(_roi_pixels == 0) / np.size(_roi_pixels) > 0.01:
        _out_of_boundaries = True

    return np.mean(_roi_pixels[_non_zero_mask]), _out_of_boundaries


def roi_fibers_density_time_point(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x,
                                  _offset_y, _offset_z, _cell_id, _direction, _time_point, _group_properties=None,
                                  _print=True, _save=True):
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
        if _print:
            print('Computing:', _experiment, _series_id, _group, _cell_id, 'roi', _time_point_roi, 'direction',
                  _direction, 'tp', _time_point, sep='\t')
        _roi_fibers_density = roi_fibers_density(_experiment, _series_id, _group, _time_point,
                                                 _time_point_roi)
        if _save:
            _time_point_fibers_densities[_time_point_roi] = _roi_fibers_density
            save.fibers_densities(_experiment, _series_id, _group, _time_point, _time_point_fibers_densities)
        return _roi_fibers_density


def roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x, _offset_y,
                               _offset_z, _cell_id, _direction, _time_points, _print=True, _save=True,
                               _out_of_borders=True):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _fibers_densities = [
        roi_fibers_density_time_point(
            _experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x, _offset_y, _offset_z, _cell_id,
            _direction, _time_point, _group_properties, _print, _save
        ) for _time_point in range(min(_time_points, len(_group_properties['time_points'])))
    ]

    if not _out_of_borders:
        return longest_fibers_densities_ascending_sequence(_fibers_densities)
    else:
        return _fibers_densities


def roi_fibers_density_by_time_pairs(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x,
                                     _offset_y, _offset_z, _direction, _time_points, _print=True, _save=True,
                                     _out_of_borders=True):
    return {
        'left_cell': roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z,
                                                _offset_x, _offset_y, _offset_z, 'left_cell', _direction, _time_points,
                                                _print, _save, _out_of_borders),
        'right_cell': roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z,
                                                 _offset_x, _offset_y, _offset_z, 'right_cell', _direction,
                                                 _time_points, _print, _save, _out_of_borders),
    }


def minimum_time_points(_experiments_tuples):
    _minimum_time_points = 10000
    for _tuple in _experiments_tuples:
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        _minimum_time_points = min(_minimum_time_points, len(_group_properties['time_points']))

    return _minimum_time_points


def longest_fibers_densities_ascending_sequence(_fibers_densities):
    _out_of_boundaries = np.array([_fibers_density[1] for _fibers_density in _fibers_densities])
    if False not in _out_of_boundaries:
        return []

    _idx_pairs = np.where(np.diff(np.hstack(([False], _out_of_boundaries == False, [False]))))[0].reshape(-1, 2)
    _start_longest_seq = _idx_pairs[np.diff(_idx_pairs, axis=1).argmax(), 0]
    _end_longest_seq_indices = np.where(_out_of_boundaries[_start_longest_seq:] == 1)

    if len(_end_longest_seq_indices[0]) > 0:
        _longest_seq = _fibers_densities[_start_longest_seq:_end_longest_seq_indices[0][0] + _start_longest_seq]
    else:
        _longest_seq = _fibers_densities[_start_longest_seq:]

    return [_fibers_density[0] for _fibers_density in _longest_seq]


def longest_same_indices_shared_in_borders_sub_array(_fibers_densities1, _fibers_densities2):
    _out_of_boundaries1 = np.array([_fibers_density[1] for _fibers_density in _fibers_densities1])
    _out_of_boundaries2 = np.array([_fibers_density[1] for _fibers_density in _fibers_densities2])

    # not-and on both
    _out_of_boundaries = np.logical_not(
        np.logical_and(np.logical_not(_out_of_boundaries1), np.logical_not(_out_of_boundaries2))
    )

    _fibers_densities1_filtered = []
    _fibers_densities2_filtered = []
    for _fibers_density1, _fibers_density2, _out_of_boundaries_value in \
            zip(_fibers_densities1, _fibers_densities2, _out_of_boundaries):
        _fibers_densities1_filtered.append((_fibers_density1[0], _out_of_boundaries_value))
        _fibers_densities2_filtered.append((_fibers_density2[0], _out_of_boundaries_value))

    return longest_fibers_densities_ascending_sequence(_fibers_densities1_filtered), \
        longest_fibers_densities_ascending_sequence(_fibers_densities2_filtered)


def remove_blacklist(_experiment, _series_id, _group, _fibers_densities):
    _blacklist = load.blacklist(_experiment, _series_id, _group)

    if len(_blacklist) > 0:
        _out_of_boundaries = np.array([_fibers_density[1] for _fibers_density in _fibers_densities])
        _out_of_boundaries[list(_blacklist.keys())] = False
        return [
            (_fibers_density, _out_of_boundaries_value)
            for _fibers_density, _out_of_boundaries_value in zip(_fibers_densities, _out_of_boundaries)
        ]
    else:
        return _fibers_densities


def coordinates_mean(_coordinates):
    _x_mean = np.mean([_coordinate[0] for _coordinate in _coordinates])
    _y_mean = np.mean([_coordinate[1] for _coordinate in _coordinates])
    _z_mean = np.mean([_coordinate[2] for _coordinate in _coordinates])

    return _x_mean, _y_mean, _z_mean


def smooth_coordinates_in_time(_coordinates, _n=5):
    _coordinates_smoothed = []
    for _n_index in range(1, _n):
        _current_coordinates = _coordinates[:_n_index]
        _coordinates_smoothed.append(coordinates_mean(_current_coordinates))

    for _index in range(_n, len(_coordinates) + 1):
        _current_coordinates = _coordinates[_index - _n:_index]
        _coordinates_smoothed.append(coordinates_mean(_current_coordinates))

    return _coordinates_smoothed
