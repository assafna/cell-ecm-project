import os
from multiprocessing.pool import Pool

import numpy as np
from PIL import Image

from libs.config_lib import CPUS_TO_USE
from libs.experiments import load, paths, filtering, config
from libs.experiments.compute import roi_by_microns
from libs.experiments.config import ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH, ROI_START_BY_AVERAGE_CELL_DIAMETER, \
    AVERAGE_CELL_DIAMETER_IN_MICRONS

OFFSET_X = 0
OFFSET_Y = 0
OFFSET_Z = 0
DIRECTION = 'inside'


def get_roi(_experiment, _series_id, _group_properties, _time_point, _cell_id, _direction):
    if ROI_START_BY_AVERAGE_CELL_DIAMETER:
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
    else:
        if _cell_id == 'cell':
            _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
                _experiment, _series_id, _group_properties['cell_id']) * 2
        else:
            _cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(
                _experiment, _series_id, _group_properties['cells_ids'][_cell_id]) * 2
    return roi_by_microns(
        _resolution_x=_group_properties['time_points'][_time_point]['resolutions']['x'],
        _resolution_y=_group_properties['time_points'][_time_point]['resolutions']['y'],
        _resolution_z=_group_properties['time_points'][_time_point]['resolutions']['z'],
        _length_x=ROI_LENGTH,
        _length_y=ROI_HEIGHT,
        _length_z=ROI_WIDTH,
        _offset_x=OFFSET_X,
        _offset_y=OFFSET_Y,
        _offset_z=OFFSET_Z,
        _cell_coordinates=_group_properties['time_points'][_time_point][_cell_id]['coordinates'],
        _cell_diameter_in_microns=_cell_diameter_in_microns,
        _direction='right' if
        (_cell_id, _direction) == ('left_cell', 'inside') or
        (_cell_id, _direction) == ('right_cell', 'outside') or
        (_cell_id, _direction) == ('cell', 'right') else 'left'
    )


def mark_cells(_image, _group_properties, _time_point, _cell_coordinates, _cell_diameter_in_microns):
    _x, _y = _cell_coordinates['x'], _cell_coordinates['y']
    _x_cell_diameter = _cell_diameter_in_microns / _group_properties['time_points'][_time_point]['resolutions']['x']
    _y_cell_diameter = _cell_diameter_in_microns / _group_properties['time_points'][_time_point]['resolutions']['y']
    _image[int(round(_y - _y_cell_diameter / 2)):int(round(_y + _y_cell_diameter / 2)), int(round(_x))] = 255
    _image[int(round(_y)), int(round(_x - _x_cell_diameter / 2)):int(round(_x + _x_cell_diameter / 2))] = 255

    return _image


def draw_borders(_image, _roi):
    _x1, _y1, _, _x2, _y2, _ = _roi
    if 0 <= _y1 < _image.shape[0] and 0 <= _x1 < _image.shape[1] and 0 <= _x2 < _image.shape[1]:
        _image[_y1, _x1:_x2] = 255
    if 0 <= _y2 < _image.shape[0] and 0 <= _x1 < _image.shape[1] and 0 <= _x2 < _image.shape[1]:
        _image[_y2, _x1:_x2] = 255
    if 0 <= _y1 < _image.shape[0] and 0 <= _y2 < _image.shape[0] and 0 <= _x1 < _image.shape[1]:
        _image[_y1:_y2, _x1] = 255
    if 0 <= _y1 < _image.shape[0] and 0 <= _y2 < _image.shape[0] and 0 <= _x2 < _image.shape[1]:
        _image[_y1:_y2, _x2] = 255

    return _image


def process_group_pairs(_experiment, _series_id, _group, _mark_cells=True, _draw_borders=True):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _group_path = paths.images(_experiment + ' - All TPs', 'Series ' + str(_series_id), _group)
    os.makedirs(_group_path, exist_ok=True)
    for _time_point in range(0, len(_group_properties['time_points']), 1):
        print(_experiment, _series_id, _group, _time_point, sep='\t')
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _left_cell_coordinates = _group_properties['time_points'][_time_point]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][_time_point]['right_cell']['coordinates']
        _left_cell_id = _group_properties['cells_ids']['left_cell']
        _right_cell_id = _group_properties['cells_ids']['right_cell']
        # _left_cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(_experiment, _series_id, _left_cell_id) * 2
        # _right_cell_diameter_in_microns = load.mean_distance_to_surface_in_microns(_experiment, _series_id, _right_cell_id) * 2
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
        _z_cell_diameter = _cell_diameter_in_microns / _group_properties['time_points'][_time_point]['resolutions']['z']
        _z_image = _time_point_image[
                   int(round(_left_cell_coordinates['z'] - _z_cell_diameter / 2)):
                   int(round(_left_cell_coordinates['z'] + _z_cell_diameter / 2))
                   ]
        _average_across_z = np.rint(np.mean(_z_image, axis=0))

        if _mark_cells:
            _average_across_z = mark_cells(_average_across_z, _group_properties, _time_point, _left_cell_coordinates, AVERAGE_CELL_DIAMETER_IN_MICRONS)
            _average_across_z = mark_cells(_average_across_z, _group_properties, _time_point, _right_cell_coordinates, AVERAGE_CELL_DIAMETER_IN_MICRONS)

        if _draw_borders:
            _left_roi = get_roi(_experiment, _series_id, _group_properties, _time_point, 'left_cell', DIRECTION)
            _right_roi = get_roi(_experiment, _series_id, _group_properties, _time_point, 'right_cell', DIRECTION)
            _average_across_z = draw_borders(_average_across_z, _left_roi)
            _average_across_z = draw_borders(_average_across_z, _right_roi)
            # _roi = (456, 13, 381, 489, 21, 417)
            # _average_across_z = draw_borders(_average_across_z, _roi)

        # _z_image_stretched = np.repeat(
        #     _average_across_z, repeats=int(round(_average_across_z.shape[1] / _average_across_z.shape[0])), axis=0
        # )
        _img = Image.fromarray(_average_across_z).convert('L')
        _img_path = os.path.join(_group_path, str(_time_point) + '.png')
        _img.save(_img_path)


def process_group_single_cells(_experiment, _series_id, _group, _mark_cells=True, _draw_borders=True):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _group_path = paths.images(_experiment + ' - First TP - Start Not Average', 'Series ' + str(_series_id), _group)
    os.makedirs(_group_path, exist_ok=True)
    # for _time_point in range(0, len(_group_properties['time_points']), 1):
    for _time_point in [0]:
        print(_experiment, _series_id, _group, _time_point, sep='\t')
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _cell_coordinates = _group_properties['time_points'][_time_point]['cell']['coordinates']
        _cell_diameter_in_microns = AVERAGE_CELL_DIAMETER_IN_MICRONS
        _z_cell_diameter = _cell_diameter_in_microns / _group_properties['time_points'][_time_point]['resolutions']['z']
        _z_image = _time_point_image[
                   int(round(_cell_coordinates['z'] - _z_cell_diameter / 2)):
                   int(round(_cell_coordinates['z'] + _z_cell_diameter / 2))
                   ]
        _average_across_z = np.rint(np.mean(_z_image, axis=0))

        if _mark_cells:
            _average_across_z = mark_cells(_average_across_z, _group_properties, _time_point, _cell_coordinates, AVERAGE_CELL_DIAMETER_IN_MICRONS)

        if _draw_borders:
            for _direction in ['left', 'right', 'up', 'down']:
                _roi = get_roi(_experiment, _series_id, _group_properties, _time_point, 'cell', _direction)
                _average_across_z = draw_borders(_average_across_z, _roi)

        # _z_image_stretched = np.repeat(
        #     _average_across_z, repeats=int(round(_average_across_z.shape[1] / _average_across_z.shape[0])), axis=0
        # )
        _img = Image.fromarray(_average_across_z).convert('L')
        _img_path = os.path.join(_group_path, str(_time_point) + '.png')
        _img.save(_img_path)


def process_experiments(_experiments, _pairs=True):
    _tuples = []
    for _experiment in _experiments:
        _tuples += load.experiment_groups_as_tuples(_experiment)
    _p = Pool(CPUS_TO_USE)
    if _pairs:
        _answers = _p.starmap(process_group_pairs, _tuples)
    else:
        _answers = _p.starmap(process_group_single_cells, _tuples)
    _p.close()


def process_all_experiments():
    # TODO: handle single cell
    process_experiments(config.PAIRS)


if __name__ == '__main__':
    # process_all_experiments()
    # process_experiments(['SN45'])
    # process_group('SN16', 2, 'cells_0_1')
    process_experiments(['Single_Cell_Ortal'], _pairs=False)
