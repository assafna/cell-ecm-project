import math

from libs.simulations import save
from libs.compute import two_lines_intersection, is_value_between_values
from libs.simulations.config import CELL_DIAMETER, EPSILON
from libs.simulations.load import elements, fibers_lengths, intersections, properties, fibers_densities


def roi_fibers_density(_simulation, _time_point, _roi):
    _elements = elements(_simulation)
    _time_point_0_fiber_lengths = fibers_lengths(_simulation, _time_point=0)
    _time_point_fiber_lengths = fibers_lengths(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point)
    _roi_lines = [
        line(_roi[0], _roi[1], _roi[0], _roi[3]),
        line(_roi[0], _roi[3], _roi[2], _roi[3]),
        line(_roi[2], _roi[1], _roi[2], _roi[3]),
        line(_roi[0], _roi[1], _roi[2], _roi[1])
    ]
    _roi_size = (_roi[2] - _roi[0]) * (_roi[3] - _roi[1])
    _roi_sum = 0
    for _element_index, _element in enumerate(_elements):
        _i1_id, _i2_id = _element
        _x1, _y1 = _intersections[_i1_id]
        _x2, _y2 = _intersections[_i2_id]
        _original_fiber_length = _time_point_0_fiber_lengths[_element_index]
        _current_fiber_length = _time_point_fiber_lengths[_element_index]
        _roi_sum += fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _roi, _roi_lines)

    return _roi_sum / _roi_size


def roi_fibers_density_by_time(_simulation, _length_x, _length_y, _offset_x, _offset_y, _cell_id, _direction, _time_points):
    _simulation_properties = properties(_simulation)
    _fibers_densities = []
    for _time_point in range(min(_time_points, len(_simulation_properties['time_points']))):
        _time_point_properties = _simulation_properties['time_points'][_time_point]
        _time_point_roi = roi(
            _length_x=_length_x,
            _length_y=_length_y,
            _offset_x=_offset_x,
            _offset_y=_offset_y,
            _cell_coordinates=_time_point_properties[_cell_id]['coordinates'],
            _cell_diameter=_time_point_properties[_cell_id]['diameter'],
            _direction='right' if
            (_cell_id, _direction) == ('left_cell', 'inside') or
            (_cell_id, _direction) == ('right_cell', 'outside') else 'left'
        )
        print(_simulation, _cell_id, 'roi', _time_point_roi, 'direction', _direction, 'tp', _time_point)
        _time_point_fibers_densities = fibers_densities(_simulation, _time_point)
        if _time_point_roi in _time_point_fibers_densities:
            _fibers_densities.append(_time_point_fibers_densities[_time_point_roi])
        else:
            _roi_fibers_density = roi_fibers_density(_simulation, _time_point, _time_point_roi)
            _time_point_fibers_densities[_time_point_roi] = _roi_fibers_density
            save.fibers_densities(_simulation, _time_point, _time_point_fibers_densities)
            _fibers_densities.append(_roi_fibers_density)

    return _fibers_densities


def pairs_roi_fibers_density_by_time(_simulation, _length_x, _length_y, _offset_x, _offset_y, _direction, _time_points):
    return {
        'left_cell': roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=_length_x,
            _length_y=_length_y,
            _offset_x=_offset_x,
            _offset_y=_offset_y,
            _cell_id='left_cell',
            _direction=_direction,
            _time_points=_time_points
        ),
        'right_cell': roi_fibers_density_by_time(
            _simulation=_simulation,
            _length_x=_length_x,
            _length_y=_length_y,
            _offset_x=_offset_x,
            _offset_y=_offset_y,
            _cell_id='right_cell',
            _direction=_direction,
            _time_points=_time_points
        )
    }


def cells_distance(_properties):
    _c1_x = _properties['time_points'][0]['left_cell']['coordinates']['x']
    _c2_x = _properties['time_points'][0]['right_cell']['coordinates']['x']

    return round(abs(_c2_x - _c1_x) / CELL_DIAMETER, 10)


def is_in_region(_x, _y, _roi):
    return _roi[0] - EPSILON <= _x <= _roi[2] + EPSILON and _roi[1] - EPSILON <= _y <= _roi[3] + EPSILON


def roi(_length_x, _length_y, _offset_x, _offset_y, _cell_coordinates, _cell_diameter, _direction):
    if _direction == 'right':
        _x1 = round(_cell_coordinates['x'] + _cell_diameter / 2 + _offset_x, 10)
    elif _direction == 'left':
        _x1 = round(_cell_coordinates['x'] - _cell_diameter / 2 - _length_x, 10)
    else:
        raise Exception('No such direction ' + _direction)
    _x2 = round(_x1 + _length_x, 10)
    _y1 = round(_cell_coordinates['y'] - _length_y / 2 + _offset_y, 10)
    _y2 = round(_y1 + _length_y, 10)

    return _x1, _y1, _x2, _y2


def line(_x1, _y1, _x2, _y2):
    return round(_y1 - _y2, 10), round(_x2 - _x1, 10), - round(_x1 * _y2 - _x2 * _y1, 10)


def fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _roi, _roi_lines):
    # check if in region
    _i1_in = is_in_region(_x1, _y1, _roi)
    _i2_in = is_in_region(_x2, _y2, _roi)

    # both in region
    if _i1_in and _i2_in:
        return _original_fiber_length
    else:
        # fiber line
        _f = line(_x1, _y1, _x2, _y2)

        # get intersections
        _is = [two_lines_intersection(_f, _roi_l) for _roi_l in _roi_lines]

        # add only not false and valid
        _intersections = [
            _i for _i in _is if _i and
                                is_value_between_values(_x1, _x2, _i[0]) and
                                is_value_between_values(_y1, _y2, _i[1]) and
                                is_in_region(_i[0], _i[1], _roi)
        ]

        # remove duplicates
        _intersections = list(set(_intersections))

        # one in region
        if len(_intersections) == 1:
            _px, _py = _intersections[0]
            _inside_length = math.hypot(_px - _x1, _py - _y1)\
                if is_in_region(_x1, _y1, _roi) else math.hypot(_px - _x2, _py - _y2)
            return _inside_length * (_original_fiber_length / _current_fiber_length)

        # both out region and valid
        elif len(_intersections) == 2:
            _inside_length = math.hypot(_intersections[1][0] - _intersections[0][0],
                                        _intersections[1][1] - _intersections[0][1])
            return _inside_length * (_original_fiber_length / _current_fiber_length)

    return 0