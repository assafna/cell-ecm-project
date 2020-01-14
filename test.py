import cv2

from libs.experiments import load, compute
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS


def main2():
    _experiment = 'SN41'
    _series_id = 3
    _group = 'cells_0_1'
    _fibers_densities = compute.roi_fibers_density_by_time_pairs(
        _experiment, _series_id, _group, CELL_DIAMETER_IN_MICRONS, CELL_DIAMETER_IN_MICRONS / 2, CELL_DIAMETER_IN_MICRONS / 2,
        CELL_DIAMETER_IN_MICRONS / 2, 0, 0, 'inside', 240
    )
    for _fiber_index in range(len(_fibers_densities['right_cell'])):
        print(_fibers_densities['left_cell'][_fiber_index], _fibers_densities['right_cell'][_fiber_index], sep='\t')


def main():
    _experiment = 'SN41'
    _series_id = 1
    _group = 'cells_0_1'
    _properties = load.group_properties(_experiment, _series_id, _group)
    _img_array = []
    for _time_point in range(len(_properties['time_points'])):
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _z_image = _time_point_image[_properties['time_points'][_time_point]['left_cell']['coordinates']['z']]
        _height, _width = _z_image.shape
        _size = (_width, _height)
        _img_array.append(_z_image)

    out = cv2.VideoWriter('project.avi', cv2.VideoWriter_fourcc(*'DIVX'), 15, _size)

    for _i in range(len(_img_array)):
        out.write(_img_array[_i])
    out.release()


if __name__ == '__main__':
    main2()
