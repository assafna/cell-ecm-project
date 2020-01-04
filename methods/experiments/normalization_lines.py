import math

from libs.experiments import paths, load, save
from libs.experiments.config import IMAGE_SIZE_PIXELS

JUMPS = 5
FACTOR = 20
MIN_LINE_DISTANCE = 200


def main():
    _borders_points = [(i, FACTOR) for i in range(FACTOR, IMAGE_SIZE_PIXELS - FACTOR, JUMPS)] +\
                      [(i, IMAGE_SIZE_PIXELS - FACTOR - 1) for i in range(FACTOR, IMAGE_SIZE_PIXELS - FACTOR, JUMPS)] +\
                      [(FACTOR, i) for i in range(FACTOR, IMAGE_SIZE_PIXELS - FACTOR, JUMPS)] +\
                      [(IMAGE_SIZE_PIXELS - FACTOR - 1, i) for i in range(FACTOR, IMAGE_SIZE_PIXELS - FACTOR, JUMPS)]

    for _experiment in paths.folders(paths.OBJECTS):
        print('Experiment', _experiment)
        for _series in paths.folders(paths.objects(_experiment)):
            print(_series)
            _objects = load.objects_time_point_file_data(_experiment, _series, _time_point='tp_1.txt')
            _line = None
            _line_max_distance = 0
            for _index1 in range(0, len(_borders_points)):
                for _index2 in range(_index1 + 1, len(_borders_points)):
                    _x1, _y1 = _borders_points[_index1]
                    _x2, _y2 = _borders_points[_index2]
                    _line_distance = math.sqrt((_y2 - _y1)**2 + (_x2 - _x1)**2)
                    if _line_distance < MIN_LINE_DISTANCE:
                        continue
                    _y_subtract = _y2 - _y1
                    _x_subtract = _x2 - _x1
                    _xy = _x2 * _y1 - _y2 * _x1
                    _min_distance = 1000
                    for _point in _objects:
                        _x0, _y0, _ = _point
                        _distance = (abs(_y_subtract * _x0 - _x_subtract * _y0 + _xy)) / _line_distance
                        _min_distance = min(_min_distance, _distance)
                    if _min_distance > _line_max_distance:
                        _line_max_distance = _min_distance
                        _line = [_x1, _y1, _x2, _y2]

            save.normalization_line(_experiment, _series, _line)


if __name__ == '__main__':
    main()
