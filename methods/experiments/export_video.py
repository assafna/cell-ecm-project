import os
from multiprocessing.pool import Pool

import numpy as np
from PIL import Image

from libs.experiments import load, paths, filtering
from libs.experiments.compute import roi_by_microns
from libs.experiments.config import CELL_DIAMETER_IN_MICRONS, ROI_LENGTH, ROI_HEIGHT, ROI_WIDTH

MINIMUM_TIME_POINTS = 20
OFFSET_X = CELL_DIAMETER_IN_MICRONS * 0
OFFSET_Y = CELL_DIAMETER_IN_MICRONS * 0.5
OFFSET_Z = 0
CELLS_DISTANCES = [8]
DIRECTION = 'inside'


def get_roi(_group_properties, _time_point, _cell_id, _direction):
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
        _direction='right' if
        (_cell_id, _direction) == ('left_cell', 'inside') or
        (_cell_id, _direction) == ('right_cell', 'outside') or
        (_cell_id, _direction) == ('cell', 'right') else 'left'
    )


def mark_cells(_image, _group_properties, _time_point, _cell_coordinates):
    _x, _y = _cell_coordinates['x'], _cell_coordinates['y']
    _x_cell_diameter = CELL_DIAMETER_IN_MICRONS / _group_properties['time_points'][_time_point]['resolutions']['x']
    _y_cell_diameter = CELL_DIAMETER_IN_MICRONS / _group_properties['time_points'][_time_point]['resolutions']['y']
    _image[_y - int(round(_y_cell_diameter / 2)):_y + int(round(_y_cell_diameter / 2)), _x] = 255
    _image[_y, _x - int(round(_x_cell_diameter / 2)):_x + int(round(_x_cell_diameter / 2))] = 255

    return _image


def draw_borders(_image, _roi):
    _x1, _y1, _, _x2, _y2, _ = _roi
    # if _y1 < _image.shape[0] and _x1 < _image.shape[1] and _x2 < _image.shape[1]:
    _image[_y1, _x1:_x2] = 255
    # if _y2 < _image.shape[0] and _x1 < _image.shape[1] and _x2 < _image.shape[1]:
    _image[_y2, _x1:_x2] = 255
    # if _y1 < _image.shape[0] and _y2 < _image.shape[0] and _x1 < _image.shape[1]:
    _image[_y1:_y2, _x1] = 255
    # if _y1 < _image.shape[0] and _y2 < _image.shape[0] and _x2 < _image.shape[1]:
    _image[_y1:_y2, _x2] = 255

    return _image


def process_group(_experiment, _series_id, _group, _mark_cells=True, _draw_borders=True):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    _group_path = paths.images(_experiment + ' TEST2', 'Series ' + str(_series_id), _group)
    os.makedirs(_group_path, exist_ok=True)
    for _time_point in range(0, len(_group_properties['time_points']), 1):
        print(_experiment, _series_id, _group, _time_point, sep='\t')
        _time_point_image = load.structured_image(_experiment, _series_id, _group, _time_point)
        _left_cell_coordinates = _group_properties['time_points'][_time_point]['left_cell']['coordinates']
        _right_cell_coordinates = _group_properties['time_points'][_time_point]['right_cell']['coordinates']

        _z_cell_diameter = CELL_DIAMETER_IN_MICRONS / _group_properties['time_points'][_time_point]['resolutions']['z']
        _z_image = _time_point_image[
                   _left_cell_coordinates['z'] - int(round(_z_cell_diameter / 2)):
                   _left_cell_coordinates['z'] + int(round(_z_cell_diameter / 2))
                   ]
        _average_across_z = np.rint(np.mean(_z_image, axis=0))

        if _mark_cells:
            _average_across_z = mark_cells(_average_across_z, _group_properties, _time_point, _left_cell_coordinates)
            _average_across_z = mark_cells(_average_across_z, _group_properties, _time_point, _right_cell_coordinates)

        if _draw_borders:
            _left_roi = get_roi(_group_properties, _time_point, 'left_cell', DIRECTION)
            _right_roi = get_roi(_group_properties, _time_point, 'right_cell', DIRECTION)
            _average_across_z = draw_borders(_average_across_z, _left_roi)
            _average_across_z = draw_borders(_average_across_z, _right_roi)

        # _z_image_stretched = np.repeat(
        #     _average_across_z, repeats=int(round(_average_across_z.shape[1] / _average_across_z.shape[0])), axis=0
        # )
        _img = Image.fromarray(_average_across_z).convert('L')
        _img_path = os.path.join(_group_path, str(_time_point) + '.png')
        _img.save(_img_path)


def process_experiments():
    _experiments = load.experiment_groups_as_tuples('SN16')
    _experiments = filtering.by_distances(_experiments, CELLS_DISTANCES)
    _experiments = filtering.by_band(_experiments)
    _experiments = filtering.by_time_points_amount(_experiments, MINIMUM_TIME_POINTS)

    _p = Pool()
    _answers = _p.starmap(process_group, _experiments)
    _p.close()


# if __name__ == '__main__':
    # process_experiments()
    # process_group('SN16', 1, 'cells_2_3')
