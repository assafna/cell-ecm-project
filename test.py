import math

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.interpolation import rotate
from tifffile import tifffile

import libs.experiments.compute
from libs import compute_lib
from libs.experiments import load
from libs.experiments.compute import angle_between_three_points, rotate_point_around_another_point
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def main():
    _experiment = 'SN16_CZI'
    _series = 'Series 1'
    _group = 'cells_3_4'
    _channel = 0
    _tp = 1
    _image = tifffile.imread('G:\My Drive\BGU\Thesis\Cell-ECM & Cell-ECM-Cell Project\Data\Experiments\Outputs\Pairs\SN16_CZI\Series 1\cells_3_4.tif')
    _image_tp_0_channel_0 = _image[_tp - 1][_channel]
    _image_tp_0_channel_0 = np.pad(_image_tp_0_channel_0, (400, 400))
    _cell_coordinates_tracked_z = load.cell_coordinates_tracked_z_group_file_data(_experiment, _series, _group + '.txt')
    _cell_1_coordinates = [_value + 400 for _value in _cell_coordinates_tracked_z[0][0][:-1]]
    _cell_2_coordinates = [_value + 400 for _value in _cell_coordinates_tracked_z[1][0][:-1]]
    if _cell_1_coordinates[0] <= _cell_2_coordinates[0]:
        _c = (_cell_1_coordinates[0] + 1, _cell_1_coordinates[1])
        _angle = angle_between_three_points(_cell_2_coordinates, _cell_1_coordinates, _c)
    else:
        _c = (_cell_2_coordinates[0] + 1, _cell_2_coordinates[1])
        _angle = angle_between_three_points(_cell_1_coordinates, _cell_2_coordinates, _c)
    rotated = rotate(_image_tp_0_channel_0, angle=_angle, reshape=False)
    image_center = tuple(np.array(_image_tp_0_channel_0.shape[1::-1]) / 2)
    _cell_1_coordinates_new = rotate_point_around_another_point(_point=_cell_1_coordinates, _angle_in_radians=math.radians(_angle), _around_point=image_center)
    _cell_2_coordinates_new = rotate_point_around_another_point(_point=_cell_2_coordinates, _angle_in_radians=math.radians(_angle), _around_point=image_center)
    _roi = libs.experiments.compute.roi_by_microns(_resolution_x=0.41513, _resolution_y=0.41513, _length_x=CELL_DIAMETER_IN_MICRONS, _length_y=CELL_DIAMETER_IN_MICRONS / 2, _offset_x=0, _offset_y=0, _cell_coordinates=_cell_1_coordinates_new, _direction='right')
    for i in range(int(_cell_1_coordinates_new[1]) - 4, int(_cell_1_coordinates_new[1]) + 4):
        for j in range(int(_cell_1_coordinates_new[0]) - 4, int(_cell_1_coordinates_new[0]) + 4):
            rotated[i][j] = 255
    for i in range(int(_cell_2_coordinates_new[1]) - 4, int(_cell_2_coordinates_new[1]) + 4):
        for j in range(int(_cell_2_coordinates_new[0]) - 4, int(_cell_2_coordinates_new[0]) + 4):
            rotated[i][j] = 255
    _values = []
    for i in range(int(_roi[0]), int(_roi[2])):
        for j in range(int(_roi[1]), int(_roi[3])):
            _values.append(rotated[j][i])
            rotated[j][i] = 255
    print(len(_values), np.average(_values))
    plt.imshow(rotated)
    plt.show()


if __name__ == '__main__':
    main()
