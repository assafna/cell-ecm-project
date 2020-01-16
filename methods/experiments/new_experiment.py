import math
import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from scipy.ndimage import rotate
from tifffile import tifffile
import matplotlib.pyplot as plt

from libs import save_lib
from libs.experiments import paths, load, compute
from libs.experiments.config import FIBERS_CHANNEL_INDEX

SHOW_PLOTS = False


def process_group(_experiment, _series_id, _cells_coordinates, _cell_1_id, _cell_2_id, _series_image_by_time_points, _resolutions, _overwrite=False):
    _time_points_data = []
    _time_point = 0
    _left_cell_id = None
    _right_cell_id = None
    while _time_point is not None:
        _time_point_image = _series_image_by_time_points[_time_point]
        _cell_1_coordinates = [int(round(_value)) for _value in _cells_coordinates[_cell_1_id][_time_point]]
        _cell_2_coordinates = [int(round(_value)) for _value in _cells_coordinates[_cell_2_id][_time_point]]

        # choose left and right cells
        if _time_point == 0:
            if _cell_1_coordinates[0] <= _cell_2_coordinates[0]:
                _left_cell_id, _right_cell_id = _cell_1_id, _cell_2_id
            else:
                _right_cell_id, _left_cell_id = _cell_1_id, _cell_2_id

        # check if exists
        _group_structured_path = paths.structured(
            _experiment, 'Series ' + str(_series_id), 'cells_' + str(_cell_1_id) + '_' + str(_cell_2_id)
        )
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        _time_point_pickle_path = paths.structured(
            _experiment=_experiment,
            _series='Series ' + str(_series_id),
            _group='cells_' + str(_cell_1_id) + '_' + str(_cell_2_id),
            _time_point=str(_time_point) + '.pkl'
        )
        if os.path.isfile(_properties_json_path) and not _overwrite and os.path.isfile(_time_point_pickle_path):
            # check for more time points
            if _time_point + 1 < len(_cells_coordinates[_cell_1_id]) and \
                    _cells_coordinates[_cell_1_id][_time_point + 1] is not None and \
                    _cells_coordinates[_cell_2_id][_time_point + 1] is not None:
                _time_point += 1
            else:
                _time_point = None
            continue

        print(_experiment, 'Series ' + str(_series_id), 'Cell 1 #:', _cell_1_id, 'Cell 2 #:', _cell_2_id, 'Time point:',
              _time_point, sep='\t')

        # set coordinates
        if _left_cell_id == _cell_1_id:
            _left_cell_coordinates, _right_cell_coordinates = _cell_1_coordinates, _cell_2_coordinates
        else:
            _right_cell_coordinates, _left_cell_coordinates = _cell_1_coordinates, _cell_2_coordinates

        # compute padding
        _helper_coordinates = (_left_cell_coordinates[0] + 1, _left_cell_coordinates[1])
        _angle = compute.angle_between_three_points(
            _right_cell_coordinates, _left_cell_coordinates, _helper_coordinates
        )
        _padding_x, _padding_y = compute.axes_padding(_2d_image_shape=_time_point_image[0].shape, _angle=_angle)
        _left_cell_coordinates[0] += _padding_x
        _left_cell_coordinates[1] += _padding_y
        _right_cell_coordinates[0] += _padding_x
        _right_cell_coordinates[1] += _padding_y

        # rotate image and change axes
        _time_point_image_rotated = np.array([rotate(_z, _angle) for _z in _time_point_image])
        _time_point_image_swapped = np.swapaxes(_time_point_image_rotated, 0, 1)

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_rotated[_left_cell_coordinates[2]])
            plt.show()
            plt.imshow(_time_point_image_rotated[_right_cell_coordinates[2]])
            plt.show()

        # update coordinates
        _image_center = compute.image_center_coordinates(_image_shape=reversed(_time_point_image_rotated[0].shape))
        _left_cell_coordinates = compute.rotate_point_around_another_point(
            _point=_left_cell_coordinates,
            _angle_in_radians=math.radians(_angle),
            _around_point=_image_center
        )
        _right_cell_coordinates = compute.rotate_point_around_another_point(
            _point=_right_cell_coordinates,
            _angle_in_radians=math.radians(_angle),
            _around_point=_image_center
        )
        _fixed_y = int(round((_left_cell_coordinates[1] + _right_cell_coordinates[1]) / 2))
        # y is now z
        _left_cell_coordinates[1] = _left_cell_coordinates[2]
        _right_cell_coordinates[1] = _right_cell_coordinates[2]
        # z is now y
        _left_cell_coordinates[2] = _fixed_y
        _right_cell_coordinates[2] = _fixed_y

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_swapped[_left_cell_coordinates[2]])
            plt.show()

        # swap resolutions
        _new_resolutions = {
            'x': _resolutions['x'],
            'y': _resolutions['z'],
            'z': _resolutions['y']
        }

        # second rotate, compute padding
        _helper_coordinates = (_left_cell_coordinates[0] + 1, _left_cell_coordinates[1])
        _angle = compute.angle_between_three_points(
            _right_cell_coordinates, _left_cell_coordinates, _helper_coordinates
        )
        _padding_x, _padding_y = compute.axes_padding(_2d_image_shape=_time_point_image_swapped[0].shape, _angle=_angle)
        _left_cell_coordinates[0] += _padding_x
        _left_cell_coordinates[1] += _padding_y
        _right_cell_coordinates[0] += _padding_x
        _right_cell_coordinates[1] += _padding_y

        # rotate image
        _time_point_image_swapped_rotated = np.array([rotate(_z, _angle) for _z in _time_point_image_swapped])

        # update coordinates
        _image_center = compute.image_center_coordinates(_image_shape=reversed(_time_point_image_swapped_rotated[0].shape))
        _left_cell_coordinates = compute.rotate_point_around_another_point(
            _point=_left_cell_coordinates,
            _angle_in_radians=math.radians(_angle),
            _around_point=_image_center
        )
        _right_cell_coordinates = compute.rotate_point_around_another_point(
            _point=_right_cell_coordinates,
            _angle_in_radians=math.radians(_angle),
            _around_point=_image_center
        )
        _fixed_y = int(round((_left_cell_coordinates[1] + _right_cell_coordinates[1]) / 2))
        _left_cell_coordinates[1] = _fixed_y
        _right_cell_coordinates[1] = _fixed_y

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_swapped_rotated[_left_cell_coordinates[2]])
            plt.show()

        # update resolutions
        _new_resolutions['x'] = (_angle / 90) * _new_resolutions['y'] + ((90 - _angle) / 90) * _new_resolutions['x']
        _new_resolutions['y'] = (_angle / 90) * _new_resolutions['x'] + ((90 - _angle) / 90) * _new_resolutions['y']

        _image_z, _image_y, _image_x = _time_point_image_swapped_rotated.shape
        if not 0 <= _left_cell_coordinates[0] < _image_x or not \
                0 <= _left_cell_coordinates[1] < _image_y or not \
                0 <= _left_cell_coordinates[2] < _image_z:
            break
        if not 0 <= _right_cell_coordinates[0] < _image_x or not \
                0 <= _right_cell_coordinates[1] < _image_y or not \
                0 <= _right_cell_coordinates[2] < _image_z:
            break

        # add to array
        _time_points_data.append({
            'left_cell': {
                'coordinates': {
                    'x': _left_cell_coordinates[0],
                    'y': _left_cell_coordinates[1],
                    'z': _left_cell_coordinates[2]
                }
            },
            'right_cell': {
                'coordinates': {
                    'x': _right_cell_coordinates[0],
                    'y': _right_cell_coordinates[1],
                    'z': _right_cell_coordinates[2]
                }
            },
            'resolutions': _new_resolutions
        })

        # save to pickle
        _time_point_pickle_path = paths.structured(
            _experiment=_experiment,
            _series='Series ' + str(_series_id),
            _group='cells_' + str(_cell_1_id) + '_' + str(_cell_2_id),
            _time_point=str(_time_point) + '.pkl'
        )
        save_lib.to_pickle(_time_point_image_swapped_rotated, _time_point_pickle_path)

        # check for more time points
        if _time_point + 1 < len(_cells_coordinates[_cell_1_id]) and \
                _cells_coordinates[_cell_1_id][_time_point + 1] is not None and \
                _cells_coordinates[_cell_2_id][_time_point + 1] is not None:
            _time_point += 1
        else:
            _time_point = None

    # save properties
    _properties_data = {
        'experiment': _experiment,
        'series_id': _series_id,
        'cells_ids': {
            'left_cell': _left_cell_id,
            'right_cell': _right_cell_id
        },
        'time_points': _time_points_data,
        'band': None
    }
    _group_structured_path = paths.structured(
        _experiment, 'Series ' + str(_series_id), 'cells_' + str(_cell_1_id) + '_' + str(_cell_2_id)
    )
    _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
    save_lib.to_json(_properties_data, _properties_json_path)


def process_series(_experiment, _series_id, _overwrite=False):
    _series_name = 'Series ' + str(_series_id)
    _series_image_path = paths.serieses(_experiment, 'series_' + str(_series_id) + '_bc.tif')
    _image_properties = load.image_properties(_experiment, 'Series ' + str(_series_id))
    _series_image = tifffile.imread(_series_image_path)
    _cells_coordinates = load.cell_coordinates_tracked_series_file_data(
        _experiment, 'series_' + str(_series_id) + '.txt'
    )
    _series_image_by_time_points = [
        np.array([_z[FIBERS_CHANNEL_INDEX] for _z in _series_image[_time_point]])
        for _time_point in range(_series_image.shape[0])
    ]
    for _cell_1_id in range(len(_cells_coordinates)):
        for _cell_2_id in range(_cell_1_id + 1, len(_cells_coordinates)):
            process_group(
                _experiment=_experiment,
                _series_id=_series_id,
                _cells_coordinates=_cells_coordinates,
                _cell_1_id=_cell_1_id,
                _cell_2_id=_cell_2_id,
                _series_image_by_time_points=_series_image_by_time_points,
                _resolutions=_image_properties['resolutions'],
                _overwrite=_overwrite
            )


def process_experiment(_experiment, _overwrite=False):
    _arguments = [
        (_experiment, int(_series.split('_')[1]), _overwrite)
        for _series in paths.image_files(paths.serieses(_experiment))
    ]
    _p = Pool()
    _p.starmap(process_series, _arguments)
    _p.close()


def process_all_experiments(_overwrite=False):
    for _experiment in paths.folders(paths.SERIESES):
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    # TODO: handle single cell experiments
    # process_all_experiments()
    # process_experiment('SN41')
    process_series('SN16_CZI', 12)
