import math

import numpy as np
from bresenham import bresenham
from scipy.ndimage import rotate
from tifffile import tifffile, TiffFile
import matplotlib.pyplot as plt

from libs import compute_lib
from libs.experiments import load, compute
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def main():
    _experiment = 'SN16_CZI'
    _series = 'Series 1'
    _group = 'cells_3_4'
    _channel = 0
    _tp = 1
    _image = tifffile.imread(
        'G:\My Drive\BGU\Thesis\Cell-ECM & Cell-ECM-Cell Project\Data\Experiments\Manipulations\SN16_CZI\Series 1\series_1_bc.tif')
    _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(_experiment, 'series_' + str(_series.split()[1]) + '.txt')
    _width = _image[0][0][0].shape[1]
    _height = _image[0][0][0].shape[0]
    _diagonal = math.sqrt(_width**2 + _height**2)
    _padding_x = int(round((_diagonal - _width) / 2))
    _padding_y = int(round((_diagonal - _height) / 2))
    _cell_3 = [
        int(round(_cell_coordinates_tracked[2][0][0])) + _padding_x,
        int(round(_cell_coordinates_tracked[2][0][1])) + _padding_y,
        int(round(_cell_coordinates_tracked[2][0][2])),
    ]
    _cell_4 = [
        int(round(_cell_coordinates_tracked[3][0][0])) + _padding_x,
        int(round(_cell_coordinates_tracked[3][0][1])) + _padding_y,
        int(round(_cell_coordinates_tracked[3][0][2])),
    ]
    if _cell_3[0] <= _cell_4[0]:
        _c = (_cell_3[0] + 1, _cell_3[1])
        _angle = compute.angle_between_three_points(_cell_4, _cell_3, _c)
    else:
        _c = (_cell_4[0] + 1, _cell_4[1])
        _angle = compute.angle_between_three_points(_cell_3, _cell_4, _c)
    _image_tp_fibres = np.array([rotate(np.pad(_z[0], pad_width=((_padding_y, _padding_y), (_padding_x, _padding_x))), _angle, reshape=False) for _z in _image[0]])

    # rotate along the z
    _new_image = np.swapaxes(_image_tp_fibres, 0, 1)

    # get new coordinates
    _image_center = ((_image[0][0][0].shape[0] + _padding_y * 2) / 2, (_image[0][0][0].shape[1] + _padding_x * 2) / 2)
    _cell_3_new_coordinates = compute.rotate_point_around_another_point(_cell_3, math.radians(_angle), _image_center)
    _cell_4_new_coordinates = compute.rotate_point_around_another_point(_cell_4, math.radians(_angle), _image_center)
    _cell_3_new_x = _cell_3_new_coordinates[0]
    _cell_3_new_y = _cell_3_new_coordinates[2]
    _cell_3_new_z = int((_cell_3_new_coordinates[1] + _cell_4_new_coordinates[1]) / 2)
    # _new_image[_cell_3_new_z][_cell_3_new_y][_cell_3_new_x] = 255
    _cell_4_new_x = _cell_4_new_coordinates[0]
    _cell_4_new_y = _cell_4_new_coordinates[2]
    _cell_4_new_z = _cell_3_new_z
    _res_x = 0.41513
    _res_y = 2.0
    _res_z = 0.41513

    # second rotate
    _width = _new_image[0].shape[1]
    _height = _new_image[0].shape[0]
    _diagonal = math.sqrt(_width ** 2 + _height ** 2)
    _padding_x = int(round((_diagonal - _width) / 2))
    _padding_y = int(round((_diagonal - _height) / 2))
    _cell_3_new_x = _cell_3_new_x + _padding_x
    _cell_3_new_y = _cell_3_new_y + _padding_y
    _cell_4_new_x = _cell_4_new_x + _padding_x
    _cell_4_new_y = _cell_4_new_y + _padding_y
    if _cell_3_new_x <= _cell_4_new_x:
        _c = (_cell_3_new_x + 1, _cell_3_new_y)
        _angle = compute.angle_between_three_points((_cell_4_new_x, _cell_4_new_y), (_cell_3_new_x, _cell_3_new_y), _c)
    else:
        _c = (_cell_4_new_x + 1, _cell_4_new_y)
        _angle = compute.angle_between_three_points((_cell_3_new_x, _cell_3_new_y), (_cell_4_new_x, _cell_4_new_y), _c)
    _rotated_image = np.array([rotate(np.pad(_z, pad_width=((_padding_y, _padding_y), (_padding_x, _padding_x))), _angle, reshape=False) for _z in _new_image])

    # get new coordinates
    _image_center = ((_new_image[0].shape[0] + _padding_y * 2) / 2, (_new_image[0].shape[1] + _padding_x * 2) / 2)
    _cell_3_new_coordinates = compute.rotate_point_around_another_point((_cell_3_new_x, _cell_3_new_y, _cell_3_new_z), math.radians(_angle), _image_center)
    _cell_4_new_coordinates = compute.rotate_point_around_another_point((_cell_4_new_x, _cell_4_new_y, _cell_4_new_z), math.radians(_angle), _image_center)
    _rotated_image[_cell_3_new_coordinates[2]][_cell_3_new_coordinates[1]][_cell_3_new_coordinates[0]] = 255
    plt.imshow(_rotated_image[_cell_3_new_coordinates[2]])
    plt.show()
    print('hi')

    _roi = compute_lib.roi_by_microns(_resolution_x=_res_x, _resolution_y=_res_y, _length_x=CELL_DIAMETER_IN_MICRONS,
                                      _length_y=CELL_DIAMETER_IN_MICRONS / 2, _offset_x=0, _offset_y=0,
                                      _cell_coordinates=[_cell_3_new_x, _cell_3_new_y], _direction='right')
    _values = []
    for i in range(int(_roi[0]), int(_roi[2])):
        for j in range(int(_roi[1]), int(_roi[3])):
            _values.append(_new_image[_cell_3_new_z][j][i])
            _new_image[_cell_3_new_z][j][i] = 255
    print(len(_values), np.mean(_values))
    # for i in range(int(_cell_4_new_y) - 4, int(_cell_4_new_y) + 4):
    #     for j in range(int(_cell_4_new_x) - 4, int(_cell_4_new_x) + 4):
    #         _new_image[_cell_4_new_z][i][j] = 255
    plt.imshow(_new_image[_cell_4_new_z])
    plt.show()
    print('hi')


if __name__ == '__main__':
    main()
