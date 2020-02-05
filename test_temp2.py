import numpy as np
import matplotlib.pyplot as plt
from libs import compute_lib
from libs.experiments import load, compute, organize
from libs.experiments.compute import roi_by_microns
from libs.experiments.config import ROI_LENGTH, ROI_WIDTH, ROI_HEIGHT, CELL_DIAMETER_IN_MICRONS
from methods.experiments.export_video import draw_borders

MINIMUM_TIME_POINTS = 245
OFFSET_X = 0
OFFSET_Y = CELL_DIAMETER_IN_MICRONS * 0
OFFSET_Z = 0
DERIVATIVE = 1
CELLS_DISTANCES = [7.0]
DIRECTION = 'inside'


if __name__ == '__main__':
    _experiment, _series_id, _group = 'SN16', 19, 'cells_1_3'
    _time_point = 24
    _fibers_densities = load.fibers_densities(_experiment, _series_id, _group, _time_point)
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _time_point_roi = roi_by_microns(
        _resolution_x=_group_properties['time_points'][_time_point]['resolutions']['x'],
        _resolution_y=_group_properties['time_points'][_time_point]['resolutions']['y'],
        _resolution_z=_group_properties['time_points'][_time_point]['resolutions']['z'],
        _length_x=ROI_LENGTH,
        _length_y=ROI_HEIGHT,
        _length_z=ROI_WIDTH,
        _offset_x=OFFSET_X,
        _offset_y=OFFSET_Y,
        _offset_z=OFFSET_Z,
        _cell_coordinates=_group_properties['time_points'][_time_point]['left_cell']['coordinates'],
        _direction='right'
    )
    _f = compute.roi_fibers_density(_experiment, _series_id, _group, _time_point, _time_point_roi)
    _x1, _y1, _z1, _x2, _y2, _z2 = _time_point_roi
    _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
    _z_image = _time_point_image[_z1:_z2]
    _average_across_z = np.rint(np.mean(_z_image, axis=0))
    _average_across_z = draw_borders(_average_across_z, _time_point_roi)
    plt.imshow(_average_across_z)
    plt.show()
    print('hi')
