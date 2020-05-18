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
        ('SN45', 1, 1, 5, -80, -15),
        ('SN45', 1, 1, 6, 0, 160),
        ('SN45', 1, 3, 4, -105, 0),
        ('SN45', 1, 3, 5, -125, 0),
        ('SN45', 1, 4, 5, 80, 0),
        ('SN45', 2, 2, 3, 0, -120),
        ('SN45', 3, 1, 2, -15, -130),
        ('SN45', 4, 0, 2, 30, 105)
    ]
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(process_follow, _arguments)
    _p.close()
