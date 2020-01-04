import math

import numpy as np

from libs.experiments import load, save, paths
from libs.experiments.config import MAX_DISTANCE_CHANGE


def main():
    for _experiment in paths.folders(paths.OBJECTS_Z):
        print('Experiment', _experiment)
        _objects_data = load.objects_z_experiment_file_data(_experiment)
        for _series in _objects_data:
            print('Series', _series)
            for _group in _objects_data[_series]:
                _cell_coordinates = [[_coordinates] for _coordinates in _objects_data[_series][_group][0]]
                for _tp in _objects_data[_series][_group][1:]:
                    for _cell_index, _cell in enumerate(_cell_coordinates):
                        if _cell[-1] is not None:
                            _previous_cell_x, _previous_cell_y, _previous_cell_z = _cell[-1]
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
                save.cell_coordinates_tracked_z(_experiment, _series, _group, _cell_coordinates)


if __name__ == '__main__':
    main()
