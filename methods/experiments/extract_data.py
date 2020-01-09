import math

from libs.experiments import load, paths, config
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def main():
    for _experiment in config.PAIRS:
        for _series in paths.folders(paths.fibers_density(_experiment)):
            _series_id = _series.split()[1]
            _series_information_file = 'series_' + str(_series_id) + '_bc.txt'
            _series_information = load.information_file_data(_experiment, _series_information_file)
            _resolutions = [_line for _line in _series_information if _line.startswith('Voxel size')][0].split()[2]
            _x_resolution, _y_resolution, _z_resolution = [float(_value) for _value in _resolutions.split('x')]
            for _group in paths.folders(paths.fibers_density(_experiment, _series)):
                _group_id = _group.split('cells_')[1]
                for _z_group in paths.folders(paths.fibers_density(_experiment, _series, _group)):
                    _z_group_id = _z_group.split('cells_')[1]
                    _time_points_amount = len(
                        paths.text_files(paths.fibers_density(_experiment, _series, _group, _z_group))
                    )
                    _cell_coordinates_tracked_file = 'series_' + str(_series_id) + '.txt'
                    _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(
                        _experiment, _cell_coordinates_tracked_file
                    )
                    _cell_1_id, _cell_2_id = _group_id.split('_')
                    _cell_1_coordinates = _cell_coordinates_tracked[int(_cell_1_id) - 1]
                    _cell_2_coordinates = _cell_coordinates_tracked[int(_cell_2_id) - 1]
                    _x1, _y1, _z1 = [float(_value) for _value in _cell_1_coordinates[0]]
                    _x2, _y2, _z2 = [float(_value) for _value in _cell_2_coordinates[0]]
                    _x1, _y1, _z1 = _x1 * _x_resolution, _y1 * _y_resolution, _z1 * _z_resolution
                    _x2, _y2, _z2 = _x2 * _x_resolution, _y2 * _y_resolution, _z2 * _z_resolution
                    _distance = math.sqrt((_x1 - _x2)**2 + (_y1 - _y2)**2 + (_z1 - _z2)**2) / CELL_DIAMETER_IN_MICRONS

                    print(_experiment, _series_id, _group_id, _z_group_id, _time_points_amount, _distance, sep='\t')


if __name__ == '__main__':
    main()
