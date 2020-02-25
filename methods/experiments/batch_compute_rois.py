from itertools import product
from multiprocessing.pool import Pool

import numpy as np

from fibers_density.experiments.master_vs_slave_heatmap import VALUES_BY_CELL_DIAMETER
from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, compute, save
from libs.experiments.config import AVERAGE_CELL_DIAMETER_IN_MICRONS, ROI_START_BY_AVERAGE_CELL_DIAMETER

EXPERIMENTS = ['SN18', 'SN44']
DIRECTIONS = ['inside']
OFFSETS_X = [0, 0.5]
OFFSETS_Y = [0]
OFFSETS_Z = [0]
ROI_LENGTHS = [1]
ROI_HEIGHTS = [1]
ROI_WIDTHS = [1]
CELLS_IDS = ['left_cell', 'right_cell']


def main(_experiment, _series_id, _group, _group_properties, _time_point):
    _fibers_densities_tp = load.fibers_densities(_experiment, _series_id, _group, _time_point)
    _rois = []
    for _offset_x, _offset_y, _offset_z, _roi_length, _roi_width, _roi_height, _cell_id, _direction in \
            product(OFFSETS_X, OFFSETS_Y, OFFSETS_Z, ROI_LENGTHS, ROI_WIDTHS, ROI_HEIGHTS, CELLS_IDS, DIRECTIONS):
        if ROI_START_BY_AVERAGE_CELL_DIAMETER:
            _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
        else:
            _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
                _experiment, _series_id, _group_properties['cells_ids'][_cell_id]) * 2
        _roi = compute.roi_by_microns(
            _resolution_x=_group_properties['time_points'][_time_point]['resolutions']['x'],
            _resolution_y=_group_properties['time_points'][_time_point]['resolutions']['y'],
            _resolution_z=_group_properties['time_points'][_time_point]['resolutions']['z'],
            _length_x=_roi_length,
            _length_y=_roi_height,
            _length_z=_roi_width,
            _offset_x=_offset_x,
            _offset_y=_offset_y,
            _offset_z=_offset_z,
            _cell_coordinates=_group_properties['time_points'][_time_point][_cell_id]['coordinates'],
            _cell_diameter_in_microns=_cell_diameter_in_microns,
            _direction='right' if
            (_cell_id, _direction) == ('left_cell', 'inside') or
            (_cell_id, _direction) == ('right_cell', 'outside') or
            (_cell_id, _direction) == ('cell', 'right') else 'left'
        )
        if _roi not in _fibers_densities_tp:
            _rois.append(_roi)

    _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
    for _roi in _rois:
        _fibers_densities_tp[_roi] = compute.roi_fibers_density(
            _experiment, _series_id, _group, _time_point, _roi, _time_point_image
        )

    if len(_rois) > 0:
        print(_experiment, _series_id, _group, _time_point)
        save.fibers_densities(_experiment, _series_id, _group, _time_point, _fibers_densities_tp)


if __name__ == '__main__':
    _arguments = []
    for _tuple in load.experiments_groups_as_tuples(EXPERIMENTS):
        _experiment, _series_id, _group = _tuple
        _group_properties = load.group_properties(_experiment, _series_id, _group)
        for _time_point in range(len(_group_properties['time_points'])):
            _arguments.append((_experiment, _series_id, _group, _group_properties, _time_point))
    _p = Pool(CPUS_TO_USE)
    _answers = _p.starmap(main, _arguments)
    _p.close()
