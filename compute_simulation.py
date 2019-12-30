import computations
import save_simulation
from load_simulation import elements, fibers_lengths, intersections, properties, fibers_densities
from computations import fiber


def roi_fibers_density(_simulation, _time_point, _roi):
    _elements = elements(_simulation)
    _time_point_0_fiber_lengths = fibers_lengths(_simulation, _time_point=0)
    _time_point_fiber_lengths = fibers_lengths(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point)
    _roi_lines = [
        computations.line(_roi[0], _roi[1], _roi[0], _roi[3]),
        computations.line(_roi[0], _roi[3], _roi[2], _roi[3]),
        computations.line(_roi[2], _roi[1], _roi[2], _roi[3]),
        computations.line(_roi[0], _roi[1], _roi[2], _roi[1])
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
        _time_point_roi = computations.roi(
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
            save_simulation.fibers_densities(_simulation, _time_point, _time_point_fibers_densities)
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
