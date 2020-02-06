import math
import os

import numpy as np

from libs.experiments import load, save, paths, config
from libs.experiments.config import MAX_DISTANCE_CHANGE


def process_series(_experiment, _series_id, _overwrite=False):
    _path = paths.cell_coordinates_tracked(_experiment, 'series_' + str(_series_id) + '.txt')
    if not _overwrite and os.path.isfile(_path):
        return

    print('Creating cell coordinates tracked for:', _experiment, 'Series ' + str(_series_id), sep='\t')
    _series_data = load.objects_series_file_data(_experiment, 'Series ' + str(_series_id))
    _cell_coordinates = [[_coordinates] for _coordinates in _series_data[0]]
    _time_point = 2
    for _tp in _series_data[1:]:
        for _cell_index, _cell in enumerate(_cell_coordinates):
            if _cell[-1] is not None:
                _previous_cell_x, _previous_cell_y, _previous_cell_z = _cell[-1]

                # fix path for SN18 experiment, where at time point 11 the Z is different, so ignore it
                if _experiment == 'SN18' and _time_point == 11:
                    _distances = [math.sqrt(
                        math.pow(_previous_cell_x - _optional_cell_x, 2) +
                        math.pow(_previous_cell_y - _optional_cell_y, 2)
                    ) for _optional_cell_x, _optional_cell_y, _optional_cell_z in _tp]
                else:
                    _distances = [math.sqrt(
                        math.pow(_previous_cell_x - _optional_cell_x, 2) +
                        math.pow(_previous_cell_y - _optional_cell_y, 2) +
                        math.pow(_previous_cell_z - _optional_cell_z, 2)
                    ) for _optional_cell_x, _optional_cell_y, _optional_cell_z in _tp]

                _min_distance_index = int(np.argmin(_distances)) if len(_distances) > 0 else None
                _cell_coordinates[_cell_index].append(
                    _tp[_min_distance_index] if _min_distance_index is not None and min(
                        _distances) < MAX_DISTANCE_CHANGE else None
                )
            else:
                _cell_coordinates[_cell_index].append(None)
        _time_point += 1

    save.cell_coordinates_tracked(_experiment, _series_id, _cell_coordinates)


def process_experiment(_experiment, _overwrite=False):
    for _series in paths.folders(paths.objects(_experiment)):
        process_series(_experiment, int(_series.split()[1]), _overwrite)


def process_all_experiments(_overwrite=False):
    for _experiment in config.experiments():
        process_experiment(_experiment, _overwrite)


if __name__ == '__main__':
    process_all_experiments()
