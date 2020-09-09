import math
import os

import numpy as np
from matplotlib import pyplot as plt
from multiprocess.pool import Pool
from scipy.ndimage import rotate
from tifffile import tifffile

from libs import save_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import paths, load, compute, config, organize, filtering
from libs.experiments.config import FIBERS_CHANNEL_INDEX

SHOW_PLOTS = False
SMOOTH_AMOUNT = 0


def process_real_fake(_experiment, _series_id, _cells_coordinates, _cell_1_id, _cell_2_id, _based_on_cell_id,
                      _fake_cell_x, _fake_cell_y, _series_image_by_time_points, _resolutions, _overwrite=False):
    _time_points_data = []
    _left_cell_id = None
    _right_cell_id = None
    _time_points_amount = len([_value for _value in _cells_coordinates[_based_on_cell_id] if _value is not None])

    # smooth coordinates
    _cells_coordinates_based_on_cell_smoothed = compute.smooth_coordinates_in_time(
        [_value for _value in _cells_coordinates[_based_on_cell_id] if _value is not None], _n=SMOOTH_AMOUNT
    )
    _cells_coordinates_fake_cell_smoothed = compute.smooth_coordinates_in_time(
        [_value for _value in _cells_coordinates[_based_on_cell_id] if _value is not None], _n=SMOOTH_AMOUNT
    )

    _group = 'cells_' + str(_cell_1_id) + '_' + str(_cell_2_id) + '_real_' + str(_based_on_cell_id) + \
             '_fake_' + str(_based_on_cell_id)

    # check if needed (missing time-point / properties file)
    if not _overwrite:
        _missing = False
        for _time_point in range(_time_points_amount):
            _time_point_pickle_path = paths.structured(
                _experiment=_experiment,
                _series='Series ' + str(_series_id),
                _group=_group,
                _time_point=str(_time_point) + '.pkl'
            )
            if not os.path.isfile(_time_point_pickle_path):
                _missing = True
                break
        _group_structured_path = paths.structured(
            _experiment, 'Series ' + str(_series_id), _group
        )
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        if not os.path.isfile(_properties_json_path):
            _missing = True
        if not _missing:
            return

    # compute change
    _x_change = _fake_cell_x - _cells_coordinates_based_on_cell_smoothed[0][0]
    _y_change = _fake_cell_y - _cells_coordinates_based_on_cell_smoothed[0][1]

    # running for each time point
    for _time_point in range(_time_points_amount):
        _time_point_image = _series_image_by_time_points[_time_point]
        _based_on_cell_coordinates = [_value for _value in _cells_coordinates_based_on_cell_smoothed[_time_point]]
        _fake_cell_coordinates = [_value for _value in _cells_coordinates_fake_cell_smoothed[_time_point]]

        # update coordinates of fake
        _fake_cell_coordinates = [
            _based_on_cell_coordinates[0] + _x_change,
            _based_on_cell_coordinates[1] + _y_change,
            _based_on_cell_coordinates[2]
        ]

        # choose left and right cells
        if _time_point == 0:
            if _based_on_cell_coordinates[0] <= _fake_cell_coordinates[0]:
                _left_cell_id, _right_cell_id = _cell_1_id, _cell_2_id
            else:
                _right_cell_id, _left_cell_id = _cell_1_id, _cell_2_id

        print(_experiment, 'Series ' + str(_series_id), 'Cell 1 #:', _cell_1_id, 'Cell 2 #:', _cell_2_id,
              'Based on cell #:', _based_on_cell_id, 'Time point:', _time_point, sep='\t')

        # set coordinates
        if _left_cell_id == _cell_1_id:
            _left_cell_coordinates, _right_cell_coordinates = _based_on_cell_coordinates, _fake_cell_coordinates
        else:
            _right_cell_coordinates, _left_cell_coordinates = _based_on_cell_coordinates, _fake_cell_coordinates

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
            plt.imshow(_time_point_image_rotated[int(round(_left_cell_coordinates[2]))])
            plt.show()
            plt.imshow(_time_point_image_rotated[int(round(_right_cell_coordinates[2]))])
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
        _fixed_y = (_left_cell_coordinates[1] + _right_cell_coordinates[1]) / 2
        # y is now z
        _left_cell_coordinates[1] = _left_cell_coordinates[2]
        _right_cell_coordinates[1] = _right_cell_coordinates[2]
        # z is now y
        _left_cell_coordinates[2] = _fixed_y
        _right_cell_coordinates[2] = _fixed_y

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_swapped[int(round(_left_cell_coordinates[2]))])
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
        _image_center = compute.image_center_coordinates(
            _image_shape=reversed(_time_point_image_swapped_rotated[0].shape))
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
        _fixed_y = (_left_cell_coordinates[1] + _right_cell_coordinates[1]) / 2
        _left_cell_coordinates[1] = _fixed_y
        _right_cell_coordinates[1] = _fixed_y

        if SHOW_PLOTS:
            if _time_point == 0 or _time_point == 50 or _time_point == 150:
                plt.imshow(_time_point_image_swapped_rotated[int(round(_left_cell_coordinates[2]))])
                plt.show()

        # update resolutions
        _angle = abs(_angle)
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
            _group=_group,
            _time_point=str(_time_point) + '.pkl'
        )
        save_lib.to_pickle(_time_point_image_swapped_rotated, _time_point_pickle_path)

    # save properties
    _properties_data = {
        'experiment': _experiment,
        'series_id': _series_id,
        'cells_ids': {
            'left_cell': _left_cell_id,
            'right_cell': _right_cell_id
        },
        'time_points': _time_points_data,
        'band': False,
        'fake': False,
        'static': False,
        'real_fake': True
    }
    _group_structured_path = paths.structured(_experiment, 'Series ' + str(_series_id), _group)
    _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
    save_lib.to_json(_properties_data, _properties_json_path)


def points_on_circumference(_x=0, _y=0, _radius=50, _num_of_points=100):
    return [(int(_x + (math.cos(2 * math.pi / _num_of_points * _point) * _radius)),
             int(_y + (math.sin(2 * math.pi / _num_of_points * _point) * _radius))
             ) for _point in range(0, _num_of_points + 1)]


def process_group(_experiment, _series_id, _cells_coordinates, _cell_1_id, _cell_2_id, _series_image_by_time_points,
                  _resolutions, _image_properties, _overwrite=False):
    # compute pair distance
    _cell_1_x, _cell_1_y, _ = _cells_coordinates[_cell_1_id][0]
    _cell_2_x, _cell_2_y, _ = _cells_coordinates[_cell_2_id][0]
    _pair_distance = math.sqrt(((_cell_1_x - _cell_2_x) ** 2 + (_cell_1_y - _cell_2_y) ** 2))

    # for each of the cells in the pair
    for _based_on_cell_id in [_cell_1_id, _cell_2_id]:

        # get all points in a distance of a pair distance around the relevant cell
        if _based_on_cell_id == _cell_1_id:
            _points = points_on_circumference(_x=_cell_1_x, _y=_cell_1_y, _radius=int(_pair_distance))
        else:
            _points = points_on_circumference(_x=_cell_2_x, _y=_cell_2_y, _radius=int(_pair_distance))

        # get all cell coordinates expect the relevant cell
        _all_cell_coordinates = []
        for _cell_id, _cell_coordinates in enumerate(_cells_coordinates):
            if _cell_id != _based_on_cell_id:
                _all_cell_coordinates.append(_cell_coordinates[0])

        # find the point with the maximal minimum distance from all cells
        _best_point = None
        _best_distance = 0
        for _point in _points:
            _point_x, _point_y = _point

            # ignore points out of the image
            _secure_distance = 25
            if _point_x < _secure_distance or _point_y < _secure_distance or \
                    _point_x > _image_properties['dimensions']['width'] - _secure_distance or \
                    _point_y > _image_properties['dimensions']['height'] - _secure_distance:
                continue

            _shortest_distance = math.inf
            for _cell_coordinates in _all_cell_coordinates:
                _cell_x, _cell_y, _ = _cell_coordinates
                _cells_distance = math.sqrt(((_point_x - _cell_x) ** 2 + (_point_y - _cell_y) ** 2))
                if _cells_distance < _shortest_distance:
                    _shortest_distance = _cells_distance

            if _shortest_distance > _best_distance:
                _best_distance = _shortest_distance
                _best_point = _point

        # in case no good point was found
        if _best_point is None:
            return

        # create fake pair
        process_real_fake(
            _experiment=_experiment,
            _series_id=_series_id,
            _cells_coordinates=_cells_coordinates,
            _cell_1_id=_cell_1_id,
            _cell_2_id=_cell_2_id,
            _based_on_cell_id=_based_on_cell_id,
            _fake_cell_x=_best_point[0],
            _fake_cell_y=_best_point[1],
            _series_image_by_time_points=_series_image_by_time_points,
            _resolutions=_image_properties['resolutions'],
            _overwrite=_overwrite
        )


def process_series(_experiment, _series_id, _overwrite=False):
    _series_image_path = paths.serieses(_experiment, 'series_' + str(_series_id) + '_bc.tif')
    _image_properties = load.image_properties(_experiment, _series_id)
    _series_image = tifffile.imread(_series_image_path)
    _cells_coordinates = load.cell_coordinates_tracked_series_file_data(
        _experiment, 'series_' + str(_series_id) + '.txt'
    )
    _series_image_by_time_points = [
        np.array([_z[FIBERS_CHANNEL_INDEX] for _z in _series_image[_time_point]])
        for _time_point in range(_series_image.shape[0])
    ]

    _tuples = load.experiment_groups_as_tuples(_experiment)
    _tuples = organize.by_experiment(_tuples)[_experiment]
    _tuples = filtering.by_real_cells(_tuples)
    _tuples = filtering.by_real_fake_cells(_tuples, _real_fake_cells=False)
    _tuples = filtering.by_series_id(_tuples, _series_id)

    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _cell_1_id, _cell_2_id = [int(_value) for _value in _group.split('_')[1:]]
        process_group(
            _experiment=_experiment,
            _series_id=_series_id,
            _cells_coordinates=_cells_coordinates,
            _cell_1_id=_cell_1_id,
            _cell_2_id=_cell_2_id,
            _series_image_by_time_points=_series_image_by_time_points,
            _resolutions=_image_properties['resolutions'],
            _image_properties=_image_properties,
            _overwrite=_overwrite
        )


def process_experiment(_experiment, _overwrite=False):
    _arguments = [
        (_experiment, int(_series.split('_')[1]), _overwrite)
        for _series in paths.image_files(paths.serieses(_experiment))
    ]
    _p = Pool(CPUS_TO_USE)
    _p.starmap(process_series, _arguments)
    _p.close()


def process_all_experiments(_overwrite=False):
    for _experiment in config.PAIRS:
        process_experiment(_experiment, _overwrite)


def main():
    process_all_experiments()


if __name__ == '__main__':
    main()
