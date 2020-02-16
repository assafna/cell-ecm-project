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
        ('SN16', 1, 0, 1, 42, 154, 18, 42, 461, 18),
        ('SN16', 1, 0, 2, 60, 484, 18, 238, 368, 18),
        ('SN16', 1, 0, 3, 33, 135, 18, 278, 64, 18),
        ('SN16', 1, 0, 4, 120, 59, 18, 326, 106, 18),
        ('SN16', 2, 0, 1, 55, 40, 18, 353, 40, 18),
        ('SN16', 2, 0, 2, 42, 46, 18, 42, 288, 18),
        ('SN16', 2, 0, 3, 108, 280, 18, 36, 478, 18),
        ('SN16', 2, 0, 4, 450, 223, 18, 450, 464, 18),
        ('SN16', 2, 0, 5, 115, 448, 18, 411, 448, 18),
        ('SN16', 3, 0, 1, 431, 66, 18, 431, 283, 18),
        ('SN16', 3, 0, 2, 321, 277, 18, 479, 472, 18),
        ('SN16', 3, 0, 3, 329, 478, 18, 453, 270, 18),
        ('SN16', 4, 0, 1, 77, 239, 18, 42, 466, 18),
        ('SN16', 4, 0, 2, 27, 274, 18, 230, 299, 18),
        ('SN16', 5, 0, 1, 53, 64, 18, 464, 64, 18),
        ('SN16', 5, 0, 2, 257, 33, 18, 32, 205, 18),
        ('SN16', 5, 0, 3, 312, 67, 18, 451, 267, 18),
        ('SN16', 5, 0, 4, 36, 373, 18, 314, 448, 18),
        ('SN16', 5, 0, 5, 470, 302, 18, 338, 486, 18),
        ('SN16', 6, 0, 1, 50, 113, 18, 50, 396, 18),
        ('SN16', 6, 0, 2, 42, 358, 18, 272, 460, 18),
        ('SN16', 6, 0, 3, 456, 42, 18, 465, 277, 18),
        ('SN16', 6, 0, 4, 48, 455, 18, 305, 455, 18),
        ('SN16', 7, 0, 1, 56, 57, 18, 56, 288, 18),
        ('SN16', 7, 0, 2, 39, 328, 18, 189, 478, 18),
        ('SN16', 7, 0, 3, 24, 485, 18, 160, 250, 18),
        ('SN16', 8, 0, 1, 447, 64, 18, 447, 291, 18),
        ('SN16', 8, 0, 2, 226, 62, 18, 453, 150, 18),
        ('SN16', 8, 0, 3, 187, 485, 18, 468, 412, 18),
        ('SN16', 8, 0, 4, 42, 183, 18, 42, 424, 18),
        ('SN16', 9, 0, 1, 321, 42, 18, 120, 140, 18),
        ('SN16', 9, 0, 2, 36, 70, 18, 232, 170, 18),
        ('SN16', 9, 0, 3, 325, 389, 18, 481, 484, 18),
        ('SN16', 10, 0, 1, 443, 54, 18, 292, 292, 18),
        ('SN16', 10, 0, 2, 121, 247, 18, 165, 486, 18),
        ('SN16', 11, 0, 1, 52, 51, 18, 338, 51, 18),
        ('SN16', 11, 0, 2, 222, 34, 18, 158, 298, 18),
        ('SN16', 11, 0, 3, 455, 195, 18, 415, 472, 18),
        ('SN16', 11, 0, 4, 256, 455, 18, 480, 374, 18),
        ('SN16', 12, 0, 1, 52, 64, 18, 52, 326, 18),
        ('SN16', 12, 0, 2, 34, 125, 18, 283, 238, 18),
        ('SN16', 12, 0, 3, 476, 169, 18, 388, 365, 18),
        ('SN16', 14, 0, 1, 62, 159, 18, 350, 64, 18),
        ('SN16', 14, 0, 2, 453, 64, 18, 453, 313, 18),
        ('SN16', 14, 0, 3, 40, 371, 18, 247, 459, 18),
        ('SN16', 14, 0, 4, 471, 330, 18, 302, 483, 18),
        ('SN16', 15, 0, 1, 79, 54, 18, 79, 326, 18),
        ('SN16', 15, 0, 2, 451, 165, 18, 451, 454, 18),
        ('SN16', 15, 0, 3, 256, 34, 18, 28, 166, 18),
        ('SN16', 16, 0, 1, 476, 287, 18, 304, 483, 18),
        ('SN16', 16, 0, 2, 26, 441, 18, 290, 489, 18),
        ('SN16', 16, 0, 3, 193, 27, 18, 101, 224, 18),
        ('SN16', 17, 0, 1, 37, 145, 18, 273, 34, 18),
        ('SN16', 17, 0, 2, 64, 456, 18, 370, 456, 18)
    ]
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(process_static, _arguments)
    _p.close()
