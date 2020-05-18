import math
import os
import sys
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from scipy.ndimage import rotate
from tqdm import tqdm

from libs import compute_lib
from libs.compute_lib import roi
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, save, paths, organize
from libs.experiments.config import AVERAGE_CELL_DIAMETER_IN_MICRONS, ROI_START_BY_AVERAGE_CELL_DIAMETER, NO_RETURN


def cells_distance_in_cell_size(_experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates):
    _image_properties = load.image_properties(_experiment, _series_id)
    _image_resolutions = _image_properties['resolutions']
    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
    _x1, _y1, _z1 = _x1 * _image_resolutions['x'], _y1 * _image_resolutions['y'], _z1 * _image_resolutions['z']
    _x2, _y2, _z2 = _x2 * _image_resolutions['x'], _y2 * _image_resolutions['y'], _z2 * _image_resolutions['z']

    return math.sqrt((_x1 - _x2) ** 2 + (_y1 - _y2) ** 2 + (_z1 - _z2) ** 2) / AVERAGE_CELL_DIAMETER_IN_MICRONS


def cells_distance_in_cell_size_time_point(_experiment, _series_id, _group, _time_point):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _left_cell_coordinates = \
        [list(_group_properties['time_points'][_time_point]['left_cell']['coordinates'].values())]
    _right_cell_coordinates = \
        [list(_group_properties['time_points'][_time_point]['right_cell']['coordinates'].values())]

    return cells_distance_in_cell_size(_experiment, _series_id, _left_cell_coordinates, _right_cell_coordinates)


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
                   _offset_z, _cell_coordinates, _cell_diameter_in_microns, _direction):
    # always by average
    _length_x_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_x) / _resolution_x
    _length_y_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_y) / _resolution_y
    _length_z_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_z) / _resolution_z
    _offset_x_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_x) / _resolution_x
    _offset_y_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_y) / _resolution_y
    _offset_z_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_z) / _resolution_z

    _x1, _y1, _x2, _y2 = [int(round(_value)) for _value in roi(
        _length_x=_length_x_in_pixels,
        _length_y=_length_y_in_pixels,
        _offset_x=_offset_x_in_pixels,
        _offset_y=_offset_y_in_pixels,
        _cell_coordinates={
            'x': _cell_coordinates['x'],
            'y': _cell_coordinates['y']
        },
        _cell_diameter=_cell_diameter_in_microns / _resolution_x,
        _direction=_direction
    )]

    _z1 = int(round(_cell_coordinates['z'] - _length_z_in_pixels / 2 + _offset_z_in_pixels, 10))
    _z2 = int(round(_z1 + _length_z_in_pixels, 10))

    return _x1, _y1, _z1, _x2, _y2, _z2


def roi_fibers_density(_experiment, _series_id, _group, _time_point, _roi, _time_point_image=None):
    if _time_point_image is None:
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
    _x1, _y1, _z1, _x2, _y2, _z2 = [int(round(_value)) for _value in _roi]

    _out_of_boundaries = False

    _z_shape, _y_shape, _x_shape = _time_point_image.shape
    if any([_x1 < 0, _y1 < 0, _z1 < 0, _x2 >= _x_shape, _y2 >= _y_shape, _z2 >= _z_shape]):
        _out_of_boundaries = True

        # fix boundaries
        _x1, _y1, _z1 = max(0, _x1), max(0, _y1), max(0, _z1)
        _x2, _y2, _z2 = min(_x2, _x_shape), min(_y2, _y_shape), min(_z2, _z_shape)

    _roi_pixels = _time_point_image[_z1:_z2, _y1:_y2, _x1:_x2]
    _non_zero_mask = np.nonzero(_roi_pixels)

    # check if more than 5% is black
    if not _out_of_boundaries and \
            (np.size(_roi_pixels) == 0 or np.count_nonzero(_roi_pixels == 0) / np.size(_roi_pixels) > 0.05):
        _out_of_boundaries = True

    # saturation
    if np.size(_roi_pixels) != 0:
        _saturation_fraction = np.count_nonzero(_roi_pixels == 255) / np.size(_roi_pixels)
    else:
        _saturation_fraction = None

    return np.mean(_roi_pixels[_non_zero_mask]), _out_of_boundaries, _saturation_fraction


def roi_fibers_density_time_point(_arguments):
    # stop if needed
    if os.path.isfile(os.path.join(paths.EXPERIMENTS, 'stop.txt')):
        return

    if 'group_properties' not in _arguments:
        _arguments['group_properties'] = \
            load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])

    _time_point_properties = _arguments['group_properties']['time_points'][_arguments['time_point']]
    _time_point_fibers_densities = load.fibers_densities(
        _arguments['experiment'], _arguments['series_id'], _arguments['group'], _arguments['time_point'])
    if ROI_START_BY_AVERAGE_CELL_DIAMETER:
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
    else:
        _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
            _experiment=_arguments['experiment'],
            _series_id=_arguments['series_id'],
            _cell_id=_arguments['group_properties']['cells_ids'][_arguments['cell_id']]) * 2 \
            if _arguments['cell_id'] != 'cell' else _arguments['group_properties']['cell_id'] * 2
    _time_point_roi = roi_by_microns(
        _resolution_x=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['x'],
        _resolution_y=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['y'],
        _resolution_z=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['z'],
        _length_x=_arguments['length_x'],
        _length_y=_arguments['length_y'],
        _length_z=_arguments['length_z'],
        _offset_x=_arguments['offset_x'],
        _offset_y=_arguments['offset_y'],
        _offset_z=_arguments['offset_z'],
        _cell_coordinates=
        _arguments['group_properties']['time_points'][_arguments['time_point']][_arguments['cell_id']]['coordinates'],
        _cell_diameter_in_microns=_cell_diameter_in_microns,
        _direction='right' if
        (_arguments['cell_id'], _arguments['direction']) == ('left_cell', 'inside') or
        (_arguments['cell_id'], _arguments['direction']) == ('right_cell', 'outside') or
        (_arguments['cell_id'], _arguments['direction']) == ('cell', 'right') else 'left'
    )
    if _time_point_roi in _time_point_fibers_densities:
        if NO_RETURN:
            return None
        else:
            return _arguments, _time_point_fibers_densities[_time_point_roi]
    else:
        if 'print' in _arguments and _arguments['print']:
            print('Computing:', _arguments['experiment'], _arguments['series_id'], _arguments['group'],
                  _arguments['cell_id'], 'roi', _time_point_roi, 'direction', _arguments['direction'], 'tp',
                  _arguments['time_point'], sep='\t')
        _roi_fibers_density = roi_fibers_density(_arguments['experiment'], _arguments['series_id'], _arguments['group'],
                                                 _arguments['time_point'], _time_point_roi)[:2]
        if 'save' in _arguments and _arguments['save']:
            _time_point_fibers_densities[_time_point_roi] = _roi_fibers_density
            save.fibers_densities(_arguments['experiment'], _arguments['series_id'], _arguments['group'],
                                  _arguments['time_point'], _time_point_fibers_densities)

        if NO_RETURN:
            return None
        else:
            return _arguments, _roi_fibers_density


def roi_fibers_density_by_time(_experiment, _series_id, _group, _length_x, _length_y, _length_z, _offset_x, _offset_y,
                               _offset_z, _cell_id, _direction, _time_points=sys.maxsize, _print=False, _save=True,
                               _out_of_borders=True):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _fibers_densities = [
        roi_fibers_density_time_point({
            'experiment': _experiment,
            'series_id': _series_id,
            'group': _group,
            'length_x': _length_x,
            'length_y': _length_y,
            'length_z': _length_z,
            'offset_x': _offset_x,
            'offset_y': _offset_y,
            'offset_z': _offset_z,
            'cell_id': _cell_id,
            'direction': _direction,
            'time_point': _time_point,
            'group_properties': _group_properties,
            'print': _print,
            'save': _save
        })[1] for _time_point in range(min(_time_points, len(_group_properties['time_points'])))
    ]

    if not _out_of_borders:
        return longest_fibers_densities_ascending_sequence(_fibers_densities)
    else:
        return _fibers_densities


def roi_fibers_density_by_time_pairs(_arguments):
    # stop if needed
    if os.path.isfile(os.path.join(paths.EXPERIMENTS, 'stop.txt')):
        return

    return _arguments, {
        'left_cell': roi_fibers_density_by_time(
            _arguments['experiment'],
            _arguments['series_id'],
            _arguments['group'],
            _arguments['length_x'],
            _arguments['length_y'],
            _arguments['length_z'],
            _arguments['offset_x'],
            _arguments['offset_y'],
            _arguments['offset_z'],
            'left_cell',
            _arguments['direction'],
            _arguments['time_points'] if 'time_points' in _arguments else sys.maxsize,
            _arguments['print'] if 'print' in _arguments else False,
            _arguments['save'] if 'save' in _arguments else True,
            _arguments['out_of_borders'] if 'out_of_borders' in _arguments else True
        ),
        'right_cell': roi_fibers_density_by_time(
            _arguments['experiment'],
            _arguments['series_id'],
            _arguments['group'],
            _arguments['length_x'],
            _arguments['length_y'],
            _arguments['length_z'],
            _arguments['offset_x'],
            _arguments['offset_y'],
            _arguments['offset_z'],
            'right_cell',
            _arguments['direction'],
            _arguments['time_points'] if 'time_points' in _arguments else sys.maxsize,
            _arguments['print'] if 'print' in _arguments else False,
            _arguments['save'] if 'save' in _arguments else True,
            _arguments['out_of_borders'] if 'out_of_borders' in _arguments else True
        )
    }


def minimum_time_points(_experiments_tuples):
    _minimum_time_points = sys.maxsize
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

    # based on the shortest
    _min_size = min(len(_out_of_boundaries1), len(_out_of_boundaries2))
    _out_of_boundaries1 = _out_of_boundaries1[:_min_size]
    _out_of_boundaries2 = _out_of_boundaries2[:_min_size]

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


def remove_blacklist(_experiment, _series_id, _cell_id, _fibers_densities):
    _blacklist = load.blacklist(_experiment, _series_id)

    if _cell_id in _blacklist or None in _blacklist:
        _out_of_boundaries = np.array([_fibers_density[1] for _fibers_density in _fibers_densities])
        if _cell_id in _blacklist:
            for _time_point in _blacklist[_cell_id]:
                if _time_point < len(_out_of_boundaries):
                    _out_of_boundaries[_time_point] = True
        if None in _blacklist:
            for _time_point in _blacklist[None]:
                if _time_point < len(_out_of_boundaries):
                    _out_of_boundaries[_time_point] = True
        return [
            (_fibers_density[0], _out_of_boundaries_value)
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
    if _n == 0:
        return _coordinates

    _coordinates_smoothed = []
    for _n_index in range(1, _n):
        _current_coordinates = _coordinates[:_n_index]
        _coordinates_smoothed.append(coordinates_mean(_current_coordinates))

    for _index in range(_n, len(_coordinates) + 1):
        _current_coordinates = _coordinates[_index - _n:_index]
        _coordinates_smoothed.append(coordinates_mean(_current_coordinates))

    return _coordinates_smoothed


def rois_fibers_densities(_tuple):
    _key, _rois, _saturation = _tuple
    _experiment, _series_id, _group, _time_point = _key
    _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
    return {
        (_experiment, _series_id, _group, _time_point, _roi):
            roi_fibers_density(_experiment, _series_id, _group, _time_point, _roi, _time_point_image) if _saturation
            else roi_fibers_density(_experiment, _series_id, _group, _time_point, _roi, _time_point_image)[:2]
        for _roi in _rois
    }


def fibers_densities(_tuples, _saturation=False):
    _organized_tuples = {}
    for _tuple in _tuples:
        _experiment, _series_id, _group, _time_point, _roi = _tuple
        _key, _value = (_experiment, _series_id, _group, _time_point), _roi
        if _key in _organized_tuples:
            _organized_tuples[_key].append(_value)
        else:
            _organized_tuples[_key] = [_value]

    _arguments = []
    for _key in _organized_tuples:
        _arguments.append((_key, _organized_tuples[_key], _saturation))

    _fibers_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _rois in tqdm(_p.imap_unordered(rois_fibers_densities, _arguments), total=len(_arguments),
                          desc='Computing Fibers Densities'):
            _fibers_densities.update(_rois)
        _p.close()
        _p.join()

    return _fibers_densities


def roi_time_point(_arguments):
    if 'group_properties' not in _arguments:
        _arguments['group_properties'] = \
            load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])
    _time_point_properties = _arguments['group_properties']['time_points'][_arguments['time_point']]
    _time_point_fibers_densities = load.fibers_densities(
        _arguments['experiment'], _arguments['series_id'], _arguments['group'], _arguments['time_point'])
    if ROI_START_BY_AVERAGE_CELL_DIAMETER:
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
    else:
        _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
            _experiment=_arguments['experiment'],
            _series_id=_arguments['series_id'],
            _cell_id=_arguments['group_properties']['cells_ids'][_arguments['cell_id']]) * 2 \
            if _arguments['cell_id'] != 'cell' else _arguments['group_properties']['cell_id'] * 2
    _time_point_roi = roi_by_microns(
        _resolution_x=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['x'],
        _resolution_y=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['y'],
        _resolution_z=_arguments['group_properties']['time_points'][_arguments['time_point']]['resolutions']['z'],
        _length_x=_arguments['length_x'],
        _length_y=_arguments['length_y'],
        _length_z=_arguments['length_z'],
        _offset_x=_arguments['offset_x'],
        _offset_y=_arguments['offset_y'],
        _offset_z=_arguments['offset_z'],
        _cell_coordinates=
        _arguments['group_properties']['time_points'][_arguments['time_point']][_arguments['cell_id']]['coordinates'],
        _cell_diameter_in_microns=_cell_diameter_in_microns,
        _direction='right' if
        (_arguments['cell_id'], _arguments['direction']) == ('left_cell', 'inside') or
        (_arguments['cell_id'], _arguments['direction']) == ('right_cell', 'outside') or
        (_arguments['cell_id'], _arguments['direction']) == ('cell', 'right') else 'left'
    )

    return _arguments['experiment'], _arguments['series_id'], _arguments['group'], _arguments['time_point'], _time_point_roi


def rois_by_time(_arguments):
    _group_properties = load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])

    # single time point
    if 'time_point' in _arguments:
        return _arguments, [roi_time_point(_arguments)]

    # multiple time points
    if 'time_points' not in _arguments:
        _arguments['time_points'] = sys.maxsize
    _rois = []
    for _time_point in range(min(_arguments['time_points'], len(_group_properties['time_points']))):
        _arguments['time_point'] = _time_point
        _rois.append(roi_time_point(_arguments))

    return _arguments, _rois


def rois(_arguments, _keys):
    _rois_to_compute = []
    _rois_dictionary = {}
    with Pool(CPUS_TO_USE) as _p:
        for _argument_keys, _value in tqdm(_p.imap_unordered(rois_by_time, _arguments), total=len(_arguments),
                                           desc='Computing Rois'):
            _key = tuple(_argument_keys[_argument] for _argument in _keys)
            _rois_dictionary[_key] = _value
            _rois_to_compute += _value
        _p.close()
        _p.join()

    return _rois_dictionary, _rois_to_compute


def cell_z_position_from_substrate(_experiment, _series_id, _cell_id, _time_point=0):
    _properties = load.image_properties(_experiment, _series_id)
    _series_z_position = _properties['position']['z']
    _cells_coordinates_tracked = \
        load.cell_coordinates_tracked_series_file_data(_experiment, 'series_' + str(_series_id) + '.txt')
    _cell_z_position = _cells_coordinates_tracked[int(_cell_id)][_time_point][2] * _properties['resolutions']['z']

    return _series_z_position + _cell_z_position


def group_mean_z_position_from_substrate(_experiment, _series_id, _group, _time_point=0):
    _, _left_cell_id, _right_cell_id = _group.split('_')
    _left_cell_z_position = cell_z_position_from_substrate(_experiment, _series_id, _left_cell_id, _time_point)
    _right_cell_z_position = cell_z_position_from_substrate(_experiment, _series_id, _right_cell_id, _time_point)

    return (_left_cell_z_position + _right_cell_z_position) / 2


