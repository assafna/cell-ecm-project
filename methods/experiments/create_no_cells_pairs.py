import numpy as np
from multiprocess.pool import Pool
from tifffile import tifffile

from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, paths
from libs.experiments.config import FIBERS_CHANNEL_INDEX
from methods.experiments.new_experiment_pairs import process_group


def process_static(_experiment, _series_id, _cell_1_id, _cell_2_id, _x1, _y1, _z1, _x2, _y2, _z2, _overwrite=False):
    _series_image_path = paths.serieses(_experiment, 'series_' + str(_series_id) + '_bc.tif')
    _image_properties = load.image_properties(_experiment, _series_id)
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


def process_follow(_experiment, _series_id, _cell_1_id, _cell_2_id, _x_change, _y_change, _z_change=0, _overwrite=False):
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
        ('SN16', 1, 0, 2, -50, 50),
        ('SN16', 2, 0, 1, 100, 100),
        ('SN16', 2, 0, 2, -50, -100),
        ('SN16', 3, 0, 1, 150, 50),
        ('SN16', 4, 1, 2, 50, -100),
        ('SN16', 4, 0, 1, -100, -50),
        ('SN16', 5, 0, 1, 50, -200),
        ('SN16', 6, 1, 2, -100, 100),
        ('SN16', 7, 0, 1, -150, 50),
        ('SN16', 8, 0, 1, -150, 100),
        ('SN16', 9, 1, 2, 100, 100),
        ('SN16', 10, 1, 2, 100, 0),
        ('SN16', 10, 1, 3, 100, -150),
        ('SN16', 11, 0, 1, -200, -100),
        ('SN16', 12, 1, 3, -150, -50),
        ('SN16', 12, 3, 4, -50, -50),
        ('SN16', 13, 0, 1, 100, 100),
        ('SN16', 14, 1, 2, 100, 200),
        ('SN16', 15, 0, 1, -200, -50),
        ('SN16', 16, 1, 2, -100, -100),
        ('SN16', 17, 1, 4, 50, 200),
        ('SN16', 17, 2, 4, 100, 200),
        ('SN16', 18, 0, 2, -150, 200),
        ('SN16', 19, 0, 1, 50, -50),
        ('SN16', 19, 0, 2, -50, 50),
        ('SN16', 19, 1, 3, 100, 50),
        ('SN16', 19, 1, 2, -50, -50),
        ('SN16', 19, 2, 3, -50, 100),
        ('SN16', 20, 3, 4, 100, -200),
        ('SN16', 20, 3, 5, 150, 0),
        ('SN16', 20, 3, 6, 150, 0),
        ('SN16', 20, 4, 5, 350, 100),
        ('SN16', 20, 4, 6, 350, 50),
        ('SN16', 21, 0, 1, -150, -150),
        ('SN16', 21, 0, 2, 50, -100),
        ('SN16', 21, 1, 2, -150, -250),
        ('SN16', 21, 1, 3, -150, -150),
        ('SN41', 1, 0, 1, -50, 50),
        ('SN41', 2, 1, 2, -20, 80),
        ('SN41', 3, 0, 1, 0, 70),
        ('SN41', 4, 0, 1, -90, 60),
        ('SN41', 5, 0, 1, 75, -15),
        ('SN41', 6, 0, 1, -100, -25),
        ('SN41', 6, 0, 2, 0, 130),
        ('SN41', 6, 1, 2, -80, 25),
        ('SN41', 7, 0, 1, 60, -50),
        ('SN41', 7, 0, 2, -80, -60),
        ('SN41', 7, 1, 2, 60, 35),
        ('SN41', 8, 1, 2, 20, 95),
        ('SN41', 9, 0, 1, 45, -60)
    ]
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(process_follow, _arguments)
    _p.close()
