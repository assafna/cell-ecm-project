import numpy as np
from multiprocess.pool import Pool
from tifffile import tifffile

from libs.experiments import load, paths
from libs.experiments.config import FIBERS_CHANNEL_INDEX
from methods.experiments.new_experiment import process_group


def process_static(_experiment, _series_id, _cell_1_id, _cell_2_id, _x1, _y1, _z1, _x2, _y2, _z2, _overwrite=False):
    _series_image_path = paths.serieses(_experiment, 'series_' + str(_series_id) + '_bc.tif')
    _image_properties = load.image_properties(_experiment, 'Series ' + str(_series_id))
    _series_image = tifffile.imread(_series_image_path)
    _cells_coordinates = [
        [(_x1, _y1, _z1) for _time_point in range(_series_image.shape[0])],
        [(_x2, _y2, _z2) for _time_point in range(_series_image.shape[0])]
    ]
    _series_image_by_time_points = [
        np.array([_z[FIBERS_CHANNEL_INDEX] for _z in _series_image[_time_point]])
        for _time_point in range(_series_image.shape[0])
    ]
    process_group(
        _experiment=_experiment,
        _series_id=_series_id,
        _cells_coordinates=_cells_coordinates,
        _cell_1_id=0,
        _cell_2_id=1,
        _series_image_by_time_points=_series_image_by_time_points,
        _resolutions=_image_properties['resolutions'],
        _real_cells=False,
        _fake_cell_1_id=_cell_1_id,
        _fake_cell_2_id=_cell_2_id,
        _overwrite=_overwrite
    )


def process_follow(_experiment, _series_id, _cell_1_id, _cell_2_id, _x_change, _y_change, _z_change, _overwrite=False):
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
    process_group(
        _experiment=_experiment,
        _series_id=_series_id,
        _cells_coordinates=_cells_coordinates,
        _cell_1_id=_cell_1_id,
        _cell_2_id=_cell_2_id,
        _series_image_by_time_points=_series_image_by_time_points,
        _resolutions=_image_properties['resolutions'],
        _real_cells=False,
        _x_change=_x_change,
        _y_change=_y_change,
        _z_change=_z_change,
        _overwrite=_overwrite
    )


if __name__ == '__main__':
    _arguments = [
        ('SN16', 1, 0, 2, 15, 68, 18, 15, 191, 18),
        ('SN16', 1, 0, 1, 15, 68, 18, 102, 17, 18),
        ('SN16', 1, 1, 2, 108, 157, 18, 15, 191, 18),
        ('SN16', 2, 0, 1, 189, 121, 18, 159, 190, 18),
        ('SN16', 2, 0, 2, 30, 191, 18, 159, 190, 18),
        ('SN16', 2, 1, 2, 30, 191, 18, 20, 102, 18),
        ('SN16', 2, 0, 3, 24, 30, 18, 153, 15, 18),
        ('SN16', 3, 0, 1, 166, 23, 18, 193, 115, 18),
        ('SN16', 3, 0, 2, 135, 195, 18, 175, 111, 18),
        ('SN16', 4, 0, 1, 41, 102, 18, 25, 196, 18),
        ('SN16', 5, 0, 1, 100, 19, 18, 197, 91, 18),
        ('SN16', 5, 0, 2, 21, 66, 18, 191, 25, 18),
        ('SN16', 5, 1, 2, 22, 39, 18, 39, 173, 18),
        ('SN16', 5, 0, 3, 23, 143, 18, 141, 189, 18)
    ]
    _p = Pool()
    _answers = _p.starmap(process_static, _arguments)
    _p.close()
