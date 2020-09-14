import math
import os

import numpy as np

from libs.experiments import load, save, paths, config

MAX_TRACKING_DISTANCE_CHANGE = 20


def process_series(_experiment, _series_id, _overwrite=False):
    _path = paths.cell_coordinates_tracked(_experiment, _series_id)
    if not _overwrite and os.path.isfile(_path):
        return

    print('Creating cell coordinates tracked for:', _experiment, 'Series ' + str(_series_id), sep='\t')
    _series_data = load.objects_series_file_data(_experiment, _series_id)
    _cell_coordinates = [[_coordinates] for _coordinates in _series_data[0]]
    _time_frame = 2
    for _time_frame in _series_data[1:]:
        for _cell_index, _cell in enumerate(_cell_coordinates):
            if _cell[-1] is not None:
                _previous_cell_x, _previous_cell_y, _previous_cell_z = _cell[-1]

                _distances = [math.sqrt(
                    (_previous_cell_x - _optional_cell_x) ** 2 +
                    (_previous_cell_y - _optional_cell_y) ** 2 +
                    (_previous_cell_z - _optional_cell_z) ** 2
                ) for _optional_cell_x, _optional_cell_y, _optional_cell_z in _time_frame]

                _min_distance_index = int(np.argmin(_distances)) if len(_distances) > 0 else None
                _cell_coordinates[_cell_index].append(
                    _time_frame[_min_distance_index] if _min_distance_index is not None and
                    min(_distances) < MAX_TRACKING_DISTANCE_CHANGE else None
                )
            else:
                _cell_coordinates[_cell_index].append(None)
        _time_frame += 1

    save.cell_coordinates_tracked(_experiment, _series_id, _cell_coordinates)


def process_experiment(_experiment, _overwrite=False):
    for _series in paths.folders(paths.objects(_experiment)):
        process_series(_experiment, int(_series.split()[1]), _overwrite)


def process_all_experiments(_overwrite=False):
    for _experiment in config.experiments():
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    process_all_experiments()
