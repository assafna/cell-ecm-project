import matplotlib.pyplot as plt
import numpy as np
from tifffile import tifffile

from libs.experiments import paths, load
from libs.experiments.config import IMAGE_FIBER_CHANNEL_INDEX


def main():
    _experiment = 'SN41'
    _series_id = 3
    _group = 'static_0_1'
    _x1, _y1, _z1 = 23, 226, 10
    _x2, _y2, _z2 = 228, 226, 10
    _time_frame = 34
    _series_image_path = paths.serieses(_experiment, _series_id)
    _image_properties = load.image_properties(_experiment, _series_id)
    _series_image = tifffile.imread(_series_image_path)
    _series_image_by_time_frames = [
        np.array([_z[IMAGE_FIBER_CHANNEL_INDEX] for _z in _series_image[_time_frame]])
        for _time_frame in range(_series_image.shape[0])
    ]
    plt.imshow(_series_image_by_time_frames[_time_frame][30])
    plt.show()


if __name__ == '__main__':
    main()
