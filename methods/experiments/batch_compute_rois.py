from itertools import product

import numpy as np

from libs.experiments import load, compute, save
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS

EXPERIMENT = 'SN16_CZI'
DIRECTION = 'inside'
OFFSETS_X = range(0, 16)
OFFSETS_Y = [0]
OFFSETS_Z = [0]
ROI_LENGTHS = [4, 8, 16]
ROI_WIDTHS = [4, 8, 16]
ROI_HEIGHTS = [4, 8, 16]
CELLS_IDS = ['left_cell', 'right_cell']


def main():
    for _tuple in load.experiment_groups_as_tuples(EXPERIMENT):
        print(_tuple)
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        for _time_point in range(len(_group_properties['time_points'])):
            _fibers_densities_tp = load.fibers_densities(_experiment, _series_id, _group, _time_point)
            _rois = []
            for _offset_x, _offset_y, _offset_z, _roi_length, _roi_width, _roi_height in \
                    product(OFFSETS_X, OFFSETS_Y, OFFSETS_Z, ROI_LENGTHS, ROI_WIDTHS, ROI_HEIGHTS):
                for _cell_id in CELLS_IDS:
                    _roi = compute.roi_by_microns(
                        _resolution_x=_group_properties['time_points'][_time_point]['resolutions']['x'],
                        _resolution_y=_group_properties['time_points'][_time_point]['resolutions']['y'],
                        _resolution_z=_group_properties['time_points'][_time_point]['resolutions']['z'],
                        _length_x=_roi_length * (CELL_DIAMETER_IN_MICRONS / 8),
                        _length_y=_roi_width * (CELL_DIAMETER_IN_MICRONS / 8),
                        _length_z=_roi_height * (CELL_DIAMETER_IN_MICRONS / 8),
                        _offset_x=_offset_x * (CELL_DIAMETER_IN_MICRONS / 8),
                        _offset_y=_offset_y * (CELL_DIAMETER_IN_MICRONS / 8),
                        _offset_z=_offset_z * (CELL_DIAMETER_IN_MICRONS / 8),
                        _cell_coordinates=_group_properties['time_points'][_time_point][_cell_id]['coordinates'],
                        _direction='right' if
                        (_cell_id, DIRECTION) == ('left_cell', 'inside') or
                        (_cell_id, DIRECTION) == ('right_cell', 'outside') or
                        (_cell_id, DIRECTION) == ('cell', 'right') else 'left'
                    )
                    if _roi not in _fibers_densities_tp:
                        _rois.append(_roi)

            _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
            for _roi in _rois:
                _x1, _y1, _z1, _x2, _y2, _z2 = _roi
                _fibers_densities_tp[_roi] = np.mean(_time_point_image[_z1:_z2, _y1:_y2, _x1:_x2])

            if len(_rois) > 0:
                save.fibers_densities(_experiment, _series_id, _group, _time_point, _fibers_densities_tp)


if __name__ == '__main__':
    main()
