import math
import os
from multiprocessing.pool import Pool

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate
from tifffile import tifffile

from libs import save_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import paths, load, compute, config
from libs.experiments.config import IMAGE_FIBER_CHANNEL_INDEX

SHOW_PLOTS = False
SMOOTH_AMOUNT = 0
CONVERT_TO_COLOR_DEPTH_8_BIT = False


# PROCESS:
# 1. Fiji script "czi_to_files_bc_concat"
# 2. Fiji script "objects_counter_all_tps"
# 3. Python script "cell_coordinates_tracked"
# 4. Python script "create_image_properties"
# 5. This Python script "new_experiment_pairs"
# 6. Python script "normalization"
# 7. Python script "export_video" for visualization
# 8. Python script "set_band_property" *** BETTER TO DO MANUALLY ***
# 9. Python script "update_blacklist"


def process_group(_experiment, _series_id, _cells_coordinates, _cell_1_id, _cell_2_id, _series_image_by_time_frames,
                  _resolutions, _real_cells=True, _fake_cell_1_id=None, _fake_cell_2_id=None,
                  _x_change=0, _y_change=0, _z_change=0, _overwrite=False):
    _time_frames_data = []
    _left_cell_id = None
    _right_cell_id = None
    _time_frames_amount = min(
        len([_value for _value in _cells_coordinates[_cell_1_id] if _value is not None]),
        len([_value for _value in _cells_coordinates[_cell_2_id] if _value is not None])
    )

    # smooth coordinates
    _cells_coordinates_cell_1_smoothed = compute.smooth_coordinates_in_time(
        [_value for _value in _cells_coordinates[_cell_1_id] if _value is not None], _n=SMOOTH_AMOUNT
    )
    _cells_coordinates_cell_2_smoothed = compute.smooth_coordinates_in_time(
        [_value for _value in _cells_coordinates[_cell_2_id] if _value is not None], _n=SMOOTH_AMOUNT
    )

    if _real_cells:
        _group = 'cells_' + str(_cell_1_id) + '_' + str(_cell_2_id)
    elif _fake_cell_1_id is None:
        _group = 'fake_' + str(_cell_1_id) + '_' + str(_cell_2_id)
    else:
        _group = 'static_' + str(_fake_cell_1_id) + '_' + str(_fake_cell_2_id)

    # check if needed (missing time-point / properties file)
    if not _overwrite:
        _missing = False
        for _time_frame in range(_time_frames_amount):
            _time_frame_pickle_path = paths.structured(_experiment, _series_id, _group, _time_frame)
            if not os.path.isfile(_time_frame_pickle_path):
                _missing = True
                break
        _group_structured_path = paths.structured(_experiment, _series_id, _group)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        if not os.path.isfile(_properties_json_path):
            _missing = True
        if not _missing:
            return

    # running for each time point
    for _time_frame in range(_time_frames_amount):
        _time_frame_image = _series_image_by_time_frames[_time_frame]
        _cell_1_coordinates = [_value for _value in _cells_coordinates_cell_1_smoothed[_time_frame]]
        _cell_2_coordinates = [_value for _value in _cells_coordinates_cell_2_smoothed[_time_frame]]

        # update coordinates if needed
        if any([_x_change != 0, _y_change != 0, _z_change != 0]):
            _cell_1_coordinates = [
                _cell_1_coordinates[0] + _x_change,
                _cell_1_coordinates[1] + _y_change,
                _cell_1_coordinates[2] + _z_change
            ]
            _cell_2_coordinates = [
                _cell_2_coordinates[0] + _x_change,
                _cell_2_coordinates[1] + _y_change,
                _cell_2_coordinates[2] + _z_change
            ]

        # choose left and right cells
        if _time_frame == 0:
            if _cell_1_coordinates[0] <= _cell_2_coordinates[0]:
                _left_cell_id, _right_cell_id = _cell_1_id, _cell_2_id
            else:
                _right_cell_id, _left_cell_id = _cell_1_id, _cell_2_id

        print(_experiment, 'Series ' + str(_series_id), 'Cell 1 #:', _cell_1_id, 'Cell 2 #:', _cell_2_id, 'Time point:',
              _time_frame, sep='\t')

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
        _padding_x, _padding_y = compute.axes_padding(_2d_image_shape=_time_frame_image[0].shape, _angle=_angle)
        _left_cell_coordinates[0] += _padding_x
        _left_cell_coordinates[1] += _padding_y
        _right_cell_coordinates[0] += _padding_x
        _right_cell_coordinates[1] += _padding_y

        # rotate image and change axes
        _time_frame_image_rotated = np.array([rotate(_z, _angle) for _z in _time_frame_image])
        _time_frame_image_swapped = np.swapaxes(_time_frame_image_rotated, 0, 1)

        if SHOW_PLOTS:
            plt.imshow(_time_frame_image_rotated[int(round(_left_cell_coordinates[2]))])
            plt.show()
            plt.imshow(_time_frame_image_rotated[int(round(_right_cell_coordinates[2]))])
            plt.show()

        # update coordinates
        _image_center = compute.image_center_coordinates(_image_shape=reversed(_time_frame_image_rotated[0].shape))
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
            plt.imshow(_time_frame_image_swapped[int(round(_left_cell_coordinates[2]))])
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
        _padding_x, _padding_y = compute.axes_padding(_2d_image_shape=_time_frame_image_swapped[0].shape, _angle=_angle)
        _left_cell_coordinates[0] += _padding_x
        _left_cell_coordinates[1] += _padding_y
        _right_cell_coordinates[0] += _padding_x
        _right_cell_coordinates[1] += _padding_y

        # rotate image
        _time_frame_image_swapped_rotated = np.array([rotate(_z, _angle) for _z in _time_frame_image_swapped])

        # convert to 8 bit color depth
        if CONVERT_TO_COLOR_DEPTH_8_BIT:
            _time_frame_image_swapped_rotated = \
                np.rint(_time_frame_image_swapped_rotated / (math.pow(2, 16) - 1) * (math.pow(2, 8) - 1)).astype(np.uint8)

        # update coordinates
        _image_center = compute.image_center_coordinates(
            _image_shape=reversed(_time_frame_image_swapped_rotated[0].shape))
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
            if _time_frame == 0 or _time_frame == 50 or _time_frame == 150:
                plt.imshow(_time_frame_image_swapped_rotated[int(round(_left_cell_coordinates[2]))])
                plt.show()

        # update resolutions
        _angle = abs(_angle)
        _new_resolution_x = (_angle / 90) * _new_resolutions['y'] + ((90 - _angle) / 90) * _new_resolutions['x']
        _new_resolution_y = (_angle / 90) * _new_resolutions['x'] + ((90 - _angle) / 90) * _new_resolutions['y']
        _new_resolutions['x'] = _new_resolution_x
        _new_resolutions['y'] = _new_resolution_y

        _image_z, _image_y, _image_x = _time_frame_image_swapped_rotated.shape
        if not 0 <= _left_cell_coordinates[0] < _image_x or not \
                0 <= _left_cell_coordinates[1] < _image_y or not \
                0 <= _left_cell_coordinates[2] < _image_z:
            break
        if not 0 <= _right_cell_coordinates[0] < _image_x or not \
                0 <= _right_cell_coordinates[1] < _image_y or not \
                0 <= _right_cell_coordinates[2] < _image_z:
            break

        # add to array
        _time_frames_data.append({
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
        _time_frame_pickle_path = paths.structured(_experiment, _series_id, _group, _time_frame)
        save_lib.to_pickle(_time_frame_image_swapped_rotated, _time_frame_pickle_path)

    # save properties
    if _real_cells:
        _band = None
        _fake = False
        _static = False
    elif _fake_cell_1_id is None:
        _based_on_properties = load.group_properties(_experiment, _series_id, 'cells_' + _group.split('fake_')[1])
        _band = _based_on_properties['band']
        _fake = True
        _static = False
    else:
        _band = False
        _fake = True
        _static = True
    _properties_data = {
        'experiment': _experiment,
        'series_id': _series_id,
        'cells_ids': {
            'left_cell': _left_cell_id,
            'right_cell': _right_cell_id
        },
        'time_points': _time_frames_data,
        'band': _band,
        'fake': _fake,
        'static': _static
    }
    _group_structured_path = paths.structured(_experiment, _series_id, _group)
    _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
    save_lib.to_json(_properties_data, _properties_json_path)


def process_series(_experiment, _series_id, _overwrite=False):
    _series_image_path = paths.serieses(_experiment, _series_id)
    _image_properties = load.image_properties(_experiment, _series_id)
    _series_image = tifffile.imread(_series_image_path)
    _cells_coordinates = load.cell_coordinates_tracked_series_file_data(_experiment, _series_id)
    _series_image_by_time_frames = [
        np.array([_z[IMAGE_FIBER_CHANNEL_INDEX] for _z in _series_image[_time_frame]])
        for _time_frame in range(_series_image.shape[0])
    ]
    for _cell_1_id in range(len(_cells_coordinates)):
        for _cell_2_id in range(_cell_1_id + 1, len(_cells_coordinates)):
            process_group(
                _experiment=_experiment,
                _series_id=_series_id,
                _cells_coordinates=_cells_coordinates,
                _cell_1_id=_cell_1_id,
                _cell_2_id=_cell_2_id,
                _series_image_by_time_frames=_series_image_by_time_frames,
                _resolutions=_image_properties['resolutions'],
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
    for _experiment in config.CELL_PAIRS:
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    process_all_experiments()
