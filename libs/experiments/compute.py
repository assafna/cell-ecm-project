import math
import sys
from multiprocessing.pool import Pool

import numpy as np
from scipy.ndimage import rotate
from tqdm import tqdm

from libs.compute_lib import window
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load
from libs.experiments.config import AVERAGE_CELL_DIAMETER_IN_MICRONS, \
    QUANTIFICATION_WINDOW_START_BY_AVERAGE_CELL_DIAMETER, \
    MAX_FRACTION_OUT_OF_BOUNDARIES_BLACK_PIXELS, MINIMUM_TIME_FRAMES_CORRELATION, DENSITY_TIME_FRAME, \
    HIGH_TEMPORAL_RESOLUTION_IN_MINUTES, QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER


def pair_distance_in_cell_size(_experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates):
    _image_properties = load.image_properties(_experiment, _series_id)
    _image_resolutions = _image_properties['resolutions']
    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
    _x1, _y1, _z1 = _x1 * _image_resolutions['x'], _y1 * _image_resolutions['y'], _z1 * _image_resolutions['z']
    _x2, _y2, _z2 = _x2 * _image_resolutions['x'], _y2 * _image_resolutions['y'], _z2 * _image_resolutions['z']

    return math.sqrt((_x1 - _x2) ** 2 + (_y1 - _y2) ** 2 + (_z1 - _z2) ** 2) / AVERAGE_CELL_DIAMETER_IN_MICRONS


def pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _left_cell_coordinates = \
        [list(_group_properties['time_points'][_time_frame]['left_cell']['coordinates'].values())]
    _right_cell_coordinates = \
        [list(_group_properties['time_points'][_time_frame]['right_cell']['coordinates'].values())]

    return pair_distance_in_cell_size(_experiment, _series_id, _left_cell_coordinates, _right_cell_coordinates)


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


def window_by_microns(_resolution_x, _resolution_y, _resolution_z, _length_x, _length_y, _length_z, _offset_x,
                      _offset_y,
                      _offset_z, _cell_coordinates, _cell_diameter_in_microns, _direction):
    # always by average
    _length_x_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_x) / _resolution_x
    _length_y_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_y) / _resolution_y
    _length_z_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _length_z) / _resolution_z
    _offset_x_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_x) / _resolution_x
    _offset_y_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_y) / _resolution_y
    _offset_z_in_pixels = (AVERAGE_CELL_DIAMETER_IN_MICRONS * _offset_z) / _resolution_z

    _x1, _y1, _x2, _y2 = [int(round(_value)) for _value in window(
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


def window_fiber_density(_experiment, _series_id, _group, _time_frame, _window, _time_frame_image=None):
    if _time_frame_image is None:
        _time_frame_image = load.structured_image(_experiment, _series_id, _group, _time_frame)
    _x1, _y1, _z1, _x2, _y2, _z2 = [int(round(_value)) for _value in _window]

    _out_of_boundaries = False

    _z_shape, _y_shape, _x_shape = _time_frame_image.shape
    if any([_x1 < 0, _y1 < 0, _z1 < 0, _x2 >= _x_shape, _y2 >= _y_shape, _z2 >= _z_shape]):
        _out_of_boundaries = True

        # fix boundaries
        _x1, _y1, _z1 = max(0, _x1), max(0, _y1), max(0, _z1)
        _x2, _y2, _z2 = min(_x2, _x_shape), min(_y2, _y_shape), min(_z2, _z_shape)

    _window_pixels = _time_frame_image[_z1:_z2, _y1:_y2, _x1:_x2]
    _non_zero_mask = np.nonzero(_window_pixels)

    # pixels to subtract
    _padding_x, _padding_y, _padding_z = 0, _y2 - _y1, _z2 - _z1
    _padding_x1, _padding_y1, _padding_z1 = max(0, _x1 - _padding_x), max(0, _y1 - _padding_y), max(0, _z1 - _padding_z)
    _padding_x2, _padding_y2, _padding_z2 = min(_x2 + _padding_x, _x_shape), min(_y2 + _padding_y, _y_shape), min(_z2 + _padding_z, _z_shape)
    _subtract_window = _time_frame_image[_padding_z1:_padding_z2, _padding_y1:_padding_y2, _padding_x1:_padding_x2].copy()
    _subtract_window[_padding_z:-_padding_z, _padding_y:-_padding_y, _padding_x:-_padding_x] = 0
    _subtract_value = np.mean(_subtract_window[np.nonzero(_subtract_window)])

    # count black pixel amount
    if not _out_of_boundaries and \
            (np.size(_window_pixels) == 0 or
             np.count_nonzero(_window_pixels == 0) / np.size(
                        _window_pixels) > MAX_FRACTION_OUT_OF_BOUNDARIES_BLACK_PIXELS):
        _out_of_boundaries = True

    # saturation amount
    if np.size(_window_pixels) != 0:
        _saturation_fraction = np.count_nonzero(_window_pixels == 255) / np.size(_window_pixels)
    else:
        _saturation_fraction = None

    # fiber density as mean of non zero pixels
    _fiber_density = np.mean(_window_pixels[_non_zero_mask])

    # subtract time frame background mean
    # _series_properties = load.image_properties(_experiment, _series_id)
    # _fiber_density -= _series_properties['time_frames'][_time_frame]['mean']

    # subtract the subtract value
    _fiber_density -= _subtract_value

    return _fiber_density, _out_of_boundaries, _saturation_fraction


def window_fiber_density_time_frame(_arguments):
    if 'group_properties' not in _arguments:
        _arguments['group_properties'] = \
            load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])

    _time_frame_properties = _arguments['group_properties']['time_points'][_arguments['time_point']]
    if QUANTIFICATION_WINDOW_START_BY_AVERAGE_CELL_DIAMETER:
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
    else:
        _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
            _experiment=_arguments['experiment'],
            _series_id=_arguments['series_id'],
            _cell_id=_arguments['group_properties']['cells_ids'][_arguments['cell_id']]) * 2 \
            if _arguments['cell_id'] != 'cell' else _arguments['group_properties']['cell_id'] * 2
    _time_frame_window = window_by_microns(
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
    if 'print' in _arguments and _arguments['print']:
        print('Computing:', _arguments['experiment'], _arguments['series_id'], _arguments['group'],
              _arguments['cell_id'], 'window', _time_frame_window, 'direction', _arguments['direction'], 'tp',
              _arguments['time_point'], sep='\t')
    _window_fiber_density = window_fiber_density(_arguments['experiment'], _arguments['series_id'], _arguments['group'],
                                                 _arguments['time_point'], _time_frame_window)[:2]

    return _arguments, _window_fiber_density


def longest_fiber_densities_ascending_sequence(_fiber_density):
    _out_of_boundaries = np.array([_fiber_density[1] for _fiber_density in _fiber_density])
    if False not in _out_of_boundaries:
        return []

    _idx_pairs = np.where(np.diff(np.hstack(([False], _out_of_boundaries == False, [False]))))[0].reshape(-1, 2)
    if len(_idx_pairs) == 0:
        return []

    _start_longest_seq = _idx_pairs[np.diff(_idx_pairs, axis=1).argmax(), 0]
    _end_longest_seq_indices = np.where(_out_of_boundaries[_start_longest_seq:] == 1)

    if len(_end_longest_seq_indices[0]) > 0:
        _longest_seq = _fiber_density[_start_longest_seq:_end_longest_seq_indices[0][0] + _start_longest_seq]
    else:
        _longest_seq = _fiber_density[_start_longest_seq:]

    return [_fiber_density[0] for _fiber_density in _longest_seq]


def longest_same_indices_shared_in_borders_sub_array(_fiber_densities_1, _fiber_densities_2):
    _out_of_boundaries_1 = np.array([_fiber_density[1] for _fiber_density in _fiber_densities_1])
    _out_of_boundaries_2 = np.array([_fiber_density[1] for _fiber_density in _fiber_densities_2])

    # based on the shortest
    _min_size = min(len(_out_of_boundaries_1), len(_out_of_boundaries_2))
    _out_of_boundaries_1 = _out_of_boundaries_1[:_min_size]
    _out_of_boundaries_2 = _out_of_boundaries_2[:_min_size]

    # not-and on both
    _out_of_boundaries = np.logical_not(
        np.logical_and(np.logical_not(_out_of_boundaries_1), np.logical_not(_out_of_boundaries_2))
    )

    _fiber_densities_1_filtered = []
    _fiber_densities_2_filtered = []
    for _fiber_density_1, _fiber_density_2, _out_of_boundaries_value in \
            zip(_fiber_densities_1, _fiber_densities_2, _out_of_boundaries):
        _fiber_densities_1_filtered.append((_fiber_density_1[0], _out_of_boundaries_value))
        _fiber_densities_2_filtered.append((_fiber_density_2[0], _out_of_boundaries_value))

    return longest_fiber_densities_ascending_sequence(_fiber_densities_1_filtered), \
        longest_fiber_densities_ascending_sequence(_fiber_densities_2_filtered)


def remove_blacklist(_experiment, _series_id, _cell_id, _fiber_densities):
    _blacklist = load.blacklist(_experiment, _series_id)

    if _cell_id in _blacklist or None in _blacklist:
        _out_of_boundaries = np.array([_fiber_density[1] for _fiber_density in _fiber_densities])
        if _cell_id in _blacklist:
            for _time_frame in _blacklist[_cell_id]:
                if _time_frame < len(_out_of_boundaries):
                    _out_of_boundaries[_time_frame] = True
        if None in _blacklist:
            for _time_frame in _blacklist[None]:
                if _time_frame < len(_out_of_boundaries):
                    _out_of_boundaries[_time_frame] = True
        return [
            (_fiber_density[0], _out_of_boundaries_value)
            for _fiber_density, _out_of_boundaries_value in zip(_fiber_densities, _out_of_boundaries)
        ]
    else:
        return _fiber_densities


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


def windows_fiber_densities(_tuple):
    _key, _windows, _saturation = _tuple
    _experiment, _series_id, _group, _time_frame = _key
    _time_frame_image = load.structured_image(_experiment, _series_id, _group, _time_frame)
    return {
        (_experiment, _series_id, _group, _time_frame, _window):
            window_fiber_density(_experiment, _series_id, _group, _time_frame, _window,
                                 _time_frame_image) if _saturation
            else window_fiber_density(_experiment, _series_id, _group, _time_frame, _window, _time_frame_image)[:2]
        for _window in _windows
    }


def fiber_densities(_tuples, _saturation=False):
    _organized_tuples = {}
    for _tuple in _tuples:
        _experiment, _series_id, _group, _time_frame, _window = _tuple
        _key, _value = (_experiment, _series_id, _group, _time_frame), _window
        if _key in _organized_tuples:
            _organized_tuples[_key].append(_value)
        else:
            _organized_tuples[_key] = [_value]

    _arguments = []
    for _key in _organized_tuples:
        _arguments.append((_key, _organized_tuples[_key], _saturation))

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _windows in tqdm(_p.imap_unordered(windows_fiber_densities, _arguments), total=len(_arguments),
                             desc='Computing fiber densities'):
            _fiber_densities.update(_windows)
        _p.close()
        _p.join()

    return _fiber_densities


def window_time_frame(_arguments):
    if 'group_properties' not in _arguments:
        _arguments['group_properties'] = \
            load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])
    _time_frame_properties = _arguments['group_properties']['time_points'][_arguments['time_point']]
    if QUANTIFICATION_WINDOW_START_BY_AVERAGE_CELL_DIAMETER:
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
    else:
        _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
            _experiment=_arguments['experiment'],
            _series_id=_arguments['series_id'],
            _cell_id=_arguments['group_properties']['cells_ids'][_arguments['cell_id']]) * 2 \
            if _arguments['cell_id'] != 'cell' else _arguments['group_properties']['cell_id'] * 2
    _time_frame_window = window_by_microns(
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

    return _arguments['experiment'], _arguments['series_id'], _arguments['group'], _arguments[
        'time_point'], _time_frame_window


def windows_by_time(_arguments):
    _group_properties = load.group_properties(_arguments['experiment'], _arguments['series_id'], _arguments['group'])

    # single time point
    if 'time_point' in _arguments:
        return _arguments, [window_time_frame(_arguments)]

    # multiple time points
    if 'time_points' not in _arguments:
        _arguments['time_points'] = sys.maxsize
    _windows = []
    for _time_frame in range(min(_arguments['time_points'], len(_group_properties['time_points']))):
        _arguments['time_point'] = _time_frame
        _windows.append(window_time_frame(_arguments))

    return _arguments, _windows


def windows(_arguments, _keys):
    _windows_to_compute = []
    _windows_dictionary = {}
    with Pool(CPUS_TO_USE) as _p:
        for _argument_keys, _value in tqdm(_p.imap_unordered(windows_by_time, _arguments), total=len(_arguments),
                                           desc='Computing windows'):
            _key = tuple(_argument_keys[_argument] for _argument in _keys)
            _windows_dictionary[_key] = _value
            _windows_to_compute += _value
        _p.close()
        _p.join()

    return _windows_dictionary, _windows_to_compute


def cell_z_position_from_substrate(_experiment, _series_id, _cell_id, _time_frame=0):
    _properties = load.image_properties(_experiment, _series_id)
    _series_z_position = _properties['position']['z']
    _cells_coordinates_tracked = \
        load.cell_coordinates_tracked_series_file_data(_experiment, _series_id)
    _cell_z_position = _cells_coordinates_tracked[int(_cell_id)][_time_frame][2] * _properties['resolutions']['z']

    return _series_z_position + _cell_z_position


def group_mean_z_position_from_substrate(_experiment, _series_id, _group, _time_frame=0):
    _, _left_cell_id, _right_cell_id = _group.split('_')
    _left_cell_z_position = cell_z_position_from_substrate(_experiment, _series_id, _left_cell_id, _time_frame)
    _right_cell_z_position = cell_z_position_from_substrate(_experiment, _series_id, _right_cell_id, _time_frame)

    return (_left_cell_z_position + _right_cell_z_position) / 2


def temporal_resolution_in_minutes(_experiment):
    _properties = load.image_properties(_experiment, _series_id=1)
    return round(_properties['frames_interval'] / 60)


def minimum_time_frames_for_correlation(_experiment):
    if temporal_resolution_in_minutes(_experiment) == HIGH_TEMPORAL_RESOLUTION_IN_MINUTES:
        return MINIMUM_TIME_FRAMES_CORRELATION['high_temporal_resolution']
    return MINIMUM_TIME_FRAMES_CORRELATION['regular_temporal_resolution']


def density_time_frame(_experiment):
    if temporal_resolution_in_minutes(_experiment) == HIGH_TEMPORAL_RESOLUTION_IN_MINUTES:
        return DENSITY_TIME_FRAME['high_temporal_resolution']
    return DENSITY_TIME_FRAME['regular_temporal_resolution']


def latest_time_frame_before_overlapping(_experiment, _series_id, _group, _offset_x):
    _properties = load.group_properties(_experiment, _series_id, _group)
    _latest_time_frame = len(_properties['time_points'])
    for _time_frame in range(len(_properties['time_points'])):
        _pair_distance = \
            pair_distance_in_cell_size_time_frame(_experiment, _series_id, _group, _time_frame)
        if _pair_distance - 1 - _offset_x * 2 < QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER * 2:
            _latest_time_frame = _time_frame - 1
            break

    return _latest_time_frame
