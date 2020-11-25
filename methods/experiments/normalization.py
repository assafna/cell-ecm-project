import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np

from libs import save_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, paths, config
from libs.experiments.config import AVERAGE_CELL_DIAMETER_IN_MICRONS

STEP_PERCENTAGE = 0.1
CELL_BORDERS_VALUES = [-2, -1, 0, 1, 2]


def is_in_window(_x1, _y1, _z1, _x2, _y2, _z2, _point_x, _point_y, _point_z):
    return all([_x1 <= _point_x <= _x2, _y1 <= _point_y <= _y2, _z1 <= _point_z <= _z2])


def process_series(_experiment, _series_id, _overwrite=False):
    _normalization_path = paths.normalization(_experiment, _series_id)
    if not _overwrite and os.path.isfile(_normalization_path):
        return

    _series_image_first_time_frame = load.series_image(_experiment, _series_id)[0]
    _image_properties = load.image_properties(_experiment, _series_id)
    _cells = load.objects_time_frame_file_data(_experiment, _series_id, _time_frame=1)
    _cell_diameter_x = AVERAGE_CELL_DIAMETER_IN_MICRONS / _image_properties['resolutions']['x']
    _cell_diameter_y = AVERAGE_CELL_DIAMETER_IN_MICRONS / _image_properties['resolutions']['y']
    _cell_diameter_z = AVERAGE_CELL_DIAMETER_IN_MICRONS / _image_properties['resolutions']['z']
    _z_shape, _y_shape, _x_shape = _series_image_first_time_frame.shape
    _z_step, _y_step, _x_step = [int(round(_value * STEP_PERCENTAGE)) for _value in [_z_shape, _y_shape, _x_shape]]

    _averages = []
    for _z, _y, _x in product(range(0, _z_shape, _z_step), range(0, _y_shape, _y_step), range(0, _x_shape, _x_step)):
        _x1 = _x - _cell_diameter_x / 2
        _x2 = _x1 + _cell_diameter_x
        _y1 = _y - _cell_diameter_y / 2
        _y2 = _y1 + _cell_diameter_y
        _z1 = _z - _cell_diameter_z / 2
        _z2 = _z1 + _cell_diameter_z
        _x1, _y1, _z1, _x2, _y2, _z2 = [int(round(_value)) for _value in [_x1, _y1, _z1, _x2, _y2, _z2]]

        # window is in borders
        if all([_x1 >= 0, _x2 < _x_shape, _y1 >= 0, _y2 < _y_shape, _z1 >= 0, _z2 < _z_shape]):

            # make sure no cells are around
            _in_window = False
            for _cell in _cells:
                _cell_x, _cell_y, _cell_z = _cell

                # cover entire box around the cell
                for _value_x, _value_y, _value_z in product(
                        CELL_BORDERS_VALUES, CELL_BORDERS_VALUES, CELL_BORDERS_VALUES
                ):
                    if is_in_window(
                            _x1, _y1, _z1, _x2, _y2, _z2,
                            _cell_x + _cell_diameter_x * _value_x,
                            _cell_y + _cell_diameter_y * _value_y,
                            _cell_z + _cell_diameter_z * _value_z
                    ):
                        _in_window = True
                        break

                # one cell is enough
                if _in_window:
                    break

            # one cell is enough
            if _in_window:
                continue

            # no cells are around, compute
            _averages.append(np.mean(_series_image_first_time_frame[_z1:_z2, _y1:_y2, _x1:_x2]))

    # compute average and std of averages
    _average, _std = np.mean(_averages), np.std(_averages)

    # save
    _normalization = {
        'average': _average,
        'std': _std
    }
    save_lib.to_json(_normalization, _normalization_path)

    print('Saved', _experiment, _series_id, sep='\t')


def process_experiment(_experiment, _overwrite=False):
    _arguments = []
    for _tuple in load.experiment_serieses_as_tuples(_experiment):
        _experiment, _series_id = _tuple
        _arguments.append((_experiment, _series_id, _overwrite))

    _p = Pool(CPUS_TO_USE)
    _p.starmap(process_series, _arguments)
    _p.close()


def process_experiments(_experiments, _overwrite=False):
    _arguments = []
    for _tuple in load.experiments_serieses_as_tuples(_experiments):
        _experiment, _series_id = _tuple
        _arguments.append((_experiment, _series_id, _overwrite))

    _p = Pool(CPUS_TO_USE)
    _p.starmap(process_series, _arguments)
    _p.close()


def process_all_experiments(_overwrite=False):
    process_experiments(config.all_experiments())


if __name__ == '__main__':
    process_all_experiments()
