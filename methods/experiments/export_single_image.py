import numpy as np
from tifffile import tifffile
import matplotlib.pyplot as plt

from libs.experiments import paths, load
from libs.experiments.config import FIBERS_CHANNEL_INDEX


def main():
    _experiment = 'SN41'
    _series_id = 3
    _group = 'static_0_1'
    _x1, _y1, _z1 = 23, 226, 10
    _x2, _y2, _z2 = 228, 226, 10
    _time_point = 34
    _series_image_path = paths.serieses(_experiment, 'series_' + str(_series_id) + '_bc.tif')
    _image_properties = load.image_properties(_experiment, _series_id)
    _series_image = tifffile.imread(_series_image_path)
    _series_image_by_time_points = [
        np.array([_z[FIBERS_CHANNEL_INDEX] for _z in _series_image[_time_point]])
        for _time_point in range(_series_image.shape[0])
    ]
    plt.imshow(_series_image_by_time_points[_time_point][30])
    plt.show()


if __name__ == '__main__':
    main()
