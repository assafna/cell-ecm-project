import math
import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage import rotate
from tifffile import tifffile
from tqdm import tqdm

from libs import save_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import paths, compute, load, config
from libs.experiments.config import FIBERS_CHANNEL_INDEX

SHOW_PLOTS = False
SMOOTH_AMOUNT = 0
DEGREES_XY = [0, 45, 90, 135]
DEGREES_Z = [0, 45, 90, 135]

# PROCESS:
# Same as in "new_experiment_pairs"


def process_group(_arguments):
    _time_points_data = []
    _time_points_amount = \
        len([_value for _value in _arguments['cell_coordinates'][_arguments['cell_id']] if _value is not None])

    # smooth coordinates
    _cells_coordinates_cell_smoothed = compute.smooth_coordinates_in_time(
        [_value for _value in _arguments['cell_coordinates'][_arguments['cell_id']] if _value is not None],
        _n=SMOOTH_AMOUNT
    )

    if _arguments['cell_type'] == 'real':
        _group = 'cell_' + str(_arguments['cell_id']) + '_' + str(_arguments['degrees_xy']) + '_' + \
                 str(_arguments['degrees_z'])
    elif _arguments['cell_type'] == 'fake':
        _group = 'fake_' + str(_arguments['cell_id']) + '_' + str(_arguments['degrees_xy']) + '_' + \
                 str(_arguments['degrees_z'])
    elif _arguments['cell_type'] == 'static':
        _group = 'static_' + str(_arguments['cell_id']) + '_' + str(_arguments['degrees_xy']) + '_' + \
                 str(_arguments['degrees_z'])
    else:
        raise Exception('No such cell type')

    # check if needed (missing time-point / properties file)
    if not _arguments['overwrite']:
        _missing = False
        for _time_point in range(_time_points_amount):
            _time_point_pickle_path = paths.structured(
                _experiment=_arguments['experiment'],
                _series='Series ' + str(_arguments['series_id']),
                _group=_group,
                _time_point=str(_time_point) + '.pkl'
            )
            if not os.path.isfile(_time_point_pickle_path):
                _missing = True
                break
        _group_structured_path = \
            paths.structured(_arguments['experiment'], 'Series ' + str(_arguments['series_id']), _group)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        if not os.path.isfile(_properties_json_path):
            _missing = True
        if not _missing:
            return

    # load image if needed
    if 'series_image_by_time_points' not in _arguments:
        _series_image_path = \
            paths.serieses(_arguments['experiment'], 'series_' + str(_arguments['series_id']) + '_bc.tif')
        _series_image = tifffile.imread(_series_image_path)
        _arguments['series_image_by_time_points'] = [
            np.array([_z[FIBERS_CHANNEL_INDEX] for _z in _series_image[_time_point]])
            for _time_point in range(_series_image.shape[0])
        ]

    # running for each time point
    for _time_point in range(_time_points_amount):
        _time_point_image = _arguments['series_image_by_time_points'][_time_point]
        _cell_coordinates = [_value for _value in _cells_coordinates_cell_smoothed[_time_point]]

        # update coordinates if needed
        if 'x_change' in _arguments:
            _cell_coordinates[0] += _arguments['x_change']
        if 'y_change' in _arguments:
            _cell_coordinates[1] += _arguments['y_change']
        if 'z_change' in _arguments:
            _cell_coordinates[2] += _arguments['z_change']

        # compute padding xy
        _padding_x, _padding_y = \
            compute.axes_padding(_2d_image_shape=_time_point_image[0].shape, _angle=_arguments['degrees_xy'])
        _cell_coordinates[0] += _padding_x
        _cell_coordinates[1] += _padding_y

        # rotate image and change axes
        _time_point_image_rotated = np.array([rotate(_z, _arguments['degrees_xy']) for _z in _time_point_image])
        _time_point_image_swapped = np.swapaxes(_time_point_image_rotated, 0, 1)

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_rotated[int(round(_cell_coordinates[2]))])
            plt.show()

        # update coordinates
        _image_center = compute.image_center_coordinates(_image_shape=reversed(_time_point_image_rotated[0].shape))
        _cell_coordinates = compute.rotate_point_around_another_point(
            _point=_cell_coordinates,
            _angle_in_radians=math.radians(_arguments['degrees_xy']),
            _around_point=_image_center
        )
        # y is now z, z is now y
        _cell_coordinates[1], _cell_coordinates[2] = _cell_coordinates[2], _cell_coordinates[1]

        if SHOW_PLOTS:
            plt.imshow(_time_point_image_swapped[int(round(_cell_coordinates[2]))])
            plt.show()

        # swap resolutions
        _new_resolutions = {
            'x': _arguments['resolutions']['x'],
            'y': _arguments['resolutions']['z'],
            'z': _arguments['resolutions']['y']
        }

        # second rotate, compute padding z
        _padding_x, _padding_y = \
            compute.axes_padding(_2d_image_shape=_time_point_image_swapped[0].shape, _angle=_arguments['degrees_z'])
        _cell_coordinates[0] += _padding_x
        _cell_coordinates[1] += _padding_y

        # rotate image
        _time_point_image_swapped_rotated = \
            np.array([rotate(_z, _arguments['degrees_z']) for _z in _time_point_image_swapped])

        # update coordinates
        _image_center = compute.image_center_coordinates(
            _image_shape=reversed(_time_point_image_swapped_rotated[0].shape))
        _cell_coordinates = compute.rotate_point_around_another_point(
            _point=_cell_coordinates,
            _angle_in_radians=math.radians(_arguments['degrees_z']),
            _around_point=_image_center
        )

        if SHOW_PLOTS:
            if _time_point == 0 or _time_point == 50 or _time_point == 150:
                plt.imshow(_time_point_image_swapped_rotated[int(round(_cell_coordinates[2]))])
                plt.show()

        # update resolutions
        _angle = abs(_arguments['degrees_z'])
        _new_resolutions['x'] = (_angle / 90) * _new_resolutions['y'] + ((90 - _angle) / 90) * _new_resolutions['x']
        _new_resolutions['y'] = (_angle / 90) * _new_resolutions['x'] + ((90 - _angle) / 90) * _new_resolutions['y']

        _image_z, _image_y, _image_x = _time_point_image_swapped_rotated.shape
        if not 0 <= _cell_coordinates[0] < _image_x or not \
                0 <= _cell_coordinates[1] < _image_y or not \
                0 <= _cell_coordinates[2] < _image_z:
            break

        # add to array
        _time_points_data.append({
            'cell': {
                'coordinates': {
                    'x': _cell_coordinates[0],
                    'y': _cell_coordinates[1],
                    'z': _cell_coordinates[2]
                }
            },
            'resolutions': _new_resolutions
        })

        # save to pickle
        _time_point_pickle_path = paths.structured(
            _experiment=_arguments['experiment'],
            _series='Series ' + str(_arguments['series_id']),
            _group=_group,
            _time_point=str(_time_point) + '.pkl'
        )
        save_lib.to_pickle(_time_point_image_swapped_rotated, _time_point_pickle_path)

    # save properties
    if _arguments['cell_type'] == 'real':
        _fake = False
        _static = False
    elif _arguments['cell_type'] == 'fake':
        _based_on_properties = \
            load.group_properties(_arguments['experiment'], _arguments['series_id'], 'cell_' + _group.split('fake_')[1])
        _fake = True
        _static = False
    elif _arguments['cell_type'] == 'static':
        _fake = True
        _static = True
    else:
        raise Exception('No such cell type')
    _properties_data = {
        'experiment': _arguments['experiment'],
        'series_id': _arguments['series_id'],
        'cell_id': _arguments['cell_id'],
        'time_points': _time_points_data,
        'fake': _fake,
        'static': _static
    }
    _group_structured_path = \
        paths.structured(_arguments['experiment'], 'Series ' + str(_arguments['series_id']), _group)
    _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
    save_lib.to_json(_properties_data, _properties_json_path)


def process_experiment(_experiment, _overwrite=False):
    _arguments = []
    for _series in paths.image_files(paths.serieses(_experiment)):
        _series_id = int(_series.split('_')[1])
        _image_properties = load.image_properties(_experiment, _series_id)
        _cells_coordinates = load.cell_coordinates_tracked_series_file_data(
            _experiment, 'series_' + str(_series_id) + '.txt'
        )
        for _cell_id in range(len(_cells_coordinates)):
            for _degrees_xy, _degrees_z in product(DEGREES_XY, DEGREES_Z):
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'cell_id': _cell_id,
                    'degrees_xy': _degrees_xy,
                    'degrees_z': _degrees_z,
                    'cell_coordinates': _cells_coordinates,
                    'cell_type': 'real',
                    'resolutions': _image_properties['resolutions'],
                    'overwrite': _overwrite
                })

    with Pool(CPUS_TO_USE) as _p:
        for _ in tqdm(_p.imap_unordered(process_group, _arguments), total=len(_arguments), desc='Creating'):
            pass
        _p.close()
        _p.join()


def process_all_experiments(_overwrite=False):
    for _experiment in config.SINGLE_CELL:
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    process_all_experiments()
