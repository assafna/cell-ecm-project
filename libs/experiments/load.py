import numpy as np
from tifffile import tifffile

from libs import load_lib
from libs.experiments import paths
from libs.experiments.config import IMAGE_FIBER_CHANNEL_INDEX
from methods.experiments.update_blacklist import add_to_blacklist


def experiment_serieses_as_tuples(_experiment):
    return [(_experiment, int(_series.split()[1])) for _series in paths.folders(paths.structured(_experiment))]


def experiments_serieses_as_tuples(_experiments):
    _tuples = []
    for _experiment in _experiments:
        _tuples += experiment_serieses_as_tuples(_experiment)

    return _tuples


def experiment_groups_as_tuples(_experiment):
    _tuples = []
    for _series in paths.folders(paths.structured(_experiment)):
        _series_id = int(_series.split()[1])
        for _group in paths.folders(paths.structured(_experiment, _series_id)):
            _tuples.append((_experiment, _series_id, _group))

    return _tuples


def experiments_groups_as_tuples(_experiments):
    _tuples = []
    for _experiment in _experiments:
        _tuples += experiment_groups_as_tuples(_experiment)

    return _tuples


def image_properties(_experiment, _series_id):
    return load_lib.from_json(paths.image_properties(_experiment, _series_id))


def group_properties(_experiment, _series_id, _group):
    return load_lib.from_json(paths.group_properties(_experiment, _series_id, _group))


def structured_image(_experiment, _series_id, _group, _time_frame):
    return load_lib.from_pickle(paths.structured(_experiment, _series_id, _group, _time_frame))


def mean_distance_to_surface_in_microns(_experiment, _series_id, _cell_id, _time_frame=0):
    _objects_file_path = paths.objects(_experiment, _series_id, _time_frame + 1)
    with open(_objects_file_path) as _f:
        _lines = _f.readlines()
        _headers = _lines[0].split('\t')
        _index = _headers.index('Mean dist. to surf. (micron)')
        return float(_lines[_cell_id + 1].split('\t')[_index])


def objects_time_frame_file_data(_experiment, _series_id, _time_frame):
    _file_path = paths.objects(_experiment, _series_id, _time_frame)
    _cells = []
    with open(_file_path) as _f:
        _lines = _f.readlines()

        # check for black image
        if _lines[0] == 'BLACK':
            for _cell_id in ['left_cell', 'right_cell']:
                add_to_blacklist(
                    _experiment=_experiment,
                    _series_id=_series_id,
                    _cell_id=_cell_id,
                    _time_frame_start=_time_frame,
                    _time_frame_end=_time_frame,
                    _reason='Black image'
                )

            # return previous time point data
            return objects_time_frame_file_data(_experiment, _series_id, _time_frame - 1)

        _headers = _lines[0].split('\t')
        _x_index, _y_index, _z_index = _headers.index('X'), _headers.index('Y'), _headers.index('Z')
        for _line in _lines[1:]:
            _line = _line.split('\t')
            _x, _y, _z = float(_line[_x_index]), float(_line[_y_index]), float(_line[_z_index])
            _cells.append((_x, _y, _z))

    return _cells


def objects_series_file_data(_experiment, _series_id):
    _objects_by_time = list([None] * len(paths.text_files(paths.objects(_experiment, _series_id))))
    for _time_frame_file in paths.text_files(paths.objects(_experiment, _series_id)):
        _time_frame = int(str(_time_frame_file.split('tp_')[1]).split('.')[0])
        _objects_by_time[_time_frame - 1] = objects_time_frame_file_data(_experiment, _series_id, _time_frame)

    return _objects_by_time


def cell_coordinates_tracked_series_file_data(_experiment, _series_id):
    _file_path = paths.cell_coordinates_tracked(_experiment, _series_id)
    _cells = []
    with open(_file_path) as _f:
        _lines = _f.readlines()
        for _line in _lines:
            _line = _line.replace('\n', '')
            _line = _line.split('\t')
            _coordinates_by_time = []
            for _coordinates in _line:
                if _coordinates == 'None':
                    _coordinates_by_time.append(None)
                else:
                    _coordinates_split = _coordinates.split()
                    _x, _y, _z = \
                        float(_coordinates_split[0]), float(_coordinates_split[1]), float(_coordinates_split[2])
                    _coordinates_by_time.append((_x, _y, _z))
            _cells.append(_coordinates_by_time)

    return _cells


def normalization_series_file_data(_experiment, _series_id):
    return load_lib.from_json(paths.normalization(_experiment, _series_id))


def series_image(_experiment, _series_id, _fiber_channel=True):
    _series_image_path = paths.serieses(_experiment, _series_id)
    _series_image = tifffile.imread(_series_image_path)

    if _fiber_channel:
        return np.array([
            np.array([_z[IMAGE_FIBER_CHANNEL_INDEX] for _z in _series_image[_time_frame]])
            for _time_frame in range(_series_image.shape[0])
        ])
    else:
        return _series_image


def blacklist(_experiment, _series_id=None):
    return load_lib.from_pickle(paths.blacklist(_experiment, _series_id))
