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
        ('SN20_Bleb_fromStart', 1, 0, 1, 280, 100),
        ('SN20_Bleb_fromStart', 1, 0, 2, 40, 150),
        ('SN20_Bleb_fromStart', 1, 1, 2, -35, 240),

        ('SN20_Bleb_fromStart', 2, 0, 1, -105, -140),

        ('SN20_Bleb_fromStart', 3, 0, 1, 365, 95),
        ('SN20_Bleb_fromStart', 3, 0, 2, 210, -245),
        ('SN20_Bleb_fromStart', 3, 1, 2, 145, -140),
        ('SN20_Bleb_fromStart', 3, 1, 3, 0, 225),

        ('SN20_Bleb_fromStart', 4, 0, 1, 130, 195),

        ('SN20_Bleb_fromStart', 5, 0, 1, 210, -160),
        ('SN20_Bleb_fromStart', 5, 0, 2, 170, 70),
        ('SN20_Bleb_fromStart', 5, 1, 2, 40, -150),
        ('SN20_Bleb_fromStart', 5, 1, 3, -110, -140),
        ('SN20_Bleb_fromStart', 5, 2, 3, 0, 160),

        ('SN20_Bleb_fromStart', 6, 0, 1, 50, -200),
        ('SN20_Bleb_fromStart', 6, 0, 2, -100, 100),
        ('SN20_Bleb_fromStart', 6, 1, 2, -80, -160),

        ('SN20_Bleb_fromStart', 7, 0, 1, -190, 60),
        ('SN20_Bleb_fromStart', 7, 0, 2, 185, 170),
        ('SN20_Bleb_fromStart', 7, 1, 2, -220, -50),

        ('SN20_Bleb_fromStart', 8, 0, 1, -150, 130),

        ('SN20_Bleb_fromStart', 9, 0, 1, -195, 140),

        ('SN20_Bleb_fromStart', 10, 0, 1, -90, -115),
        ('SN20_Bleb_fromStart', 10, 0, 2, 70, 320),
        ('SN20_Bleb_fromStart', 10, 1, 2, -40, 160),

        ('SN20_Bleb_fromStart', 11, 0, 2, -30, 305),
        ('SN20_Bleb_fromStart', 11, 1, 2, -220, -145),

        ('SN20_Bleb_fromStart', 12, 0, 1, -200, -130),
        ('SN20_Bleb_fromStart', 12, 0, 2, 30, 195),
        ('SN20_Bleb_fromStart', 12, 1, 2, -190, 40),
        ('SN20_Bleb_fromStart', 12, 1, 3, -140, 110),
        ('SN20_Bleb_fromStart', 12, 2, 3, 0, 120),

        ('SN20_Bleb_fromStart', 13, 0, 1, 85, 335),
        ('SN20_Bleb_fromStart', 13, 0, 2, -275, 265),
        ('SN20_Bleb_fromStart', 13, 0, 3, -65, 130),
        ('SN20_Bleb_fromStart', 13, 1, 2, -110, 260),
        ('SN20_Bleb_fromStart', 13, 1, 3, -170, -100),
        ('SN20_Bleb_fromStart', 13, 2, 3, 135, 20),

        ('SN20_Bleb_fromStart', 14, 0, 1, -235, 30),
        ('SN20_Bleb_fromStart', 14, 0, 2, 120, 230),
        ('SN20_Bleb_fromStart', 14, 0, 3, -230, 105),
        ('SN20_Bleb_fromStart', 14, 0, 4, 205, 35),
        ('SN20_Bleb_fromStart', 14, 1, 2, 110, -180),
        ('SN20_Bleb_fromStart', 14, 1, 3, -220, 25),
        ('SN20_Bleb_fromStart', 14, 1, 4, -150, 0),
        ('SN20_Bleb_fromStart', 14, 2, 3, 160, -130),
        ('SN20_Bleb_fromStart', 14, 2, 4, -75, 210),
        ('SN20_Bleb_fromStart', 14, 3, 4, 220, 105),

        ('SN20_Bleb_fromStart', 15, 0, 1, 0, 235),

        ('SN20_Bleb_fromStart', 16, 0, 1, 0, -225),
        ('SN20_Bleb_fromStart', 16, 0, 2, -80, 130),
        ('SN20_Bleb_fromStart', 16, 1, 2, -60, -120),

        ('SN20_Bleb_fromStart', 17, 0, 2, -180, 0),
        ('SN20_Bleb_fromStart', 17, 0, 3, 155, 0),
        ('SN20_Bleb_fromStart', 17, 1, 2, -225, -115),
        ('SN20_Bleb_fromStart', 17, 1, 3, -135, 20),

        ('SN20_Bleb_fromStart', 18, 0, 1, -110, -175),

        ('SN20_Bleb_fromStart', 19, 0, 1, 70, -150),
        ('SN20_Bleb_fromStart', 19, 1, 2, -100, 115),
        ('SN20_Bleb_fromStart', 19, 1, 3, 60, -170),
        ('SN20_Bleb_fromStart', 19, 2, 3, 135, 185),

        ('SN20_Bleb_fromStart', 20, 0, 1, 175, 20),
        ('SN20_Bleb_fromStart', 20, 0, 2, 205, -60),
        ('SN20_Bleb_fromStart', 20, 1, 2, -135, 80),
    ]
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(process_follow, _arguments)
    _p.close()
