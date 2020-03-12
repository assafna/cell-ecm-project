import math

from libs.simulations import save
from libs.compute_lib import two_lines_intersection, is_value_between_values, roi
from libs.simulations.config import CELL_DIAMETER, EPSILON, NET_DIMENSIONS
from libs.simulations.load import elements, fibers_lengths, intersections, properties, fibers_densities


def roi_fibers_density(_simulation, _time_point, _roi):
    _elements = elements(_simulation)
    _time_point_0_fibers_lengths = fibers_lengths(_simulation, _time_point=0)
    _time_point_fibers_lengths = fibers_lengths(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point)

    _master_roi = get_master_roi([_roi])
    _elements, _time_point_0_fibers_lengths, _time_point_fibers_lengths = elements_and_intersections_by_rois(
        _elements,
        _intersections,
        _time_point_0_fibers_lengths,
        _time_point_fibers_lengths,
        _master_roi
    )

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
        _original_fiber_length = _time_point_0_fibers_lengths[_element_index]
        _current_fiber_length = _time_point_fibers_lengths[_element_index]
        _fiber = fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _roi, _roi_lines)
        if _fiber is not None:
            _roi_sum += _fiber

    return _roi_sum / _roi_size


def rois_fibers_densities(_simulation, _time_point, _rois):
    _elements = elements(_simulation)
    _time_point_0_fibers_lengths = fibers_lengths(_simulation, _time_point=0)
    _time_point_fibers_lengths = fibers_lengths(_simulation, _time_point=_time_point)
    _time_point_fibers_densities = fibers_densities(_simulation, _time_point=_time_point)
    _intersections = intersections(_simulation, _time_point=_time_point)

    _rois_sums = {}
    _rois_lines = {}
    _rois_sizes = {}
    _rois_to_compute = []
    for _roi in _rois:
        if _roi in _time_point_fibers_densities:
            _rois_sums[_roi] = _time_point_fibers_densities[_roi]
        else:
            _rois_sums[_roi] = 0
            _rois_to_compute.append(_roi)
            _rois_lines[_roi] = [
                line(_roi[0], _roi[1], _roi[0], _roi[3]),
                line(_roi[0], _roi[3], _roi[2], _roi[3]),
                line(_roi[2], _roi[1], _roi[2], _roi[3]),
                line(_roi[0], _roi[1], _roi[2], _roi[1])
            ]
            _rois_sizes[_roi] = (_roi[2] - _roi[0]) * (_roi[3] - _roi[1])

    _master_roi = get_master_roi(_rois)
    _elements, _time_point_0_fibers_lengths, _time_point_fibers_lengths = elements_and_intersections_by_rois(
        _elements,
        _intersections,
        _time_point_0_fibers_lengths,
        _time_point_fibers_lengths,
        _master_roi
    )

    _elements_len = str(len(_elements))
    _rois_to_compute_len = str(len(_rois_to_compute))
    for _element_index, _element in enumerate(_elements):
        _i1_id, _i2_id = _element
        _x1, _y1 = _intersections[_i1_id]
        _x2, _y2 = _intersections[_i2_id]
        _original_fiber_length = _time_point_0_fibers_lengths[_element_index]
        _current_fiber_length = _time_point_fibers_lengths[_element_index]
        for _roi_index, _roi in enumerate(_rois_to_compute):
            print('\rElement ' + str(_element_index + 1) + '/' + _elements_len + ' ROI ' + str(_roi_index + 1) + '/' + _rois_to_compute_len, end='')
            _roi_lines = _rois_lines[_roi]
            _fiber = fiber(_x1, _y1, _x2, _y2, _original_fiber_length, _current_fiber_length, _roi, _roi_lines)
            if _fiber is not None:
                _rois_sums[_roi] += _fiber

    for _roi in _rois_to_compute:
        _rois_sums[_roi] = _rois_sums[_roi] / _rois_sizes[_roi]

    return _rois_sums


def roi_fibers_density_time_point(_arguments):
    if 'properties' not in _arguments:
        _arguments['properties'] = properties(_arguments['simulation'])
    _time_point_properties = _arguments['properties']['time_points'][_arguments['time_point']]
    _time_point_fibers_densities = fibers_densities(_arguments['simulation'], _arguments['time_point'])
    if _arguments['direction'] == 'inside':
        _new_direction = 'right' if _arguments['cell_id'] == 'left_cell' else 'left'
    elif _arguments['direction'] == 'outside':
        _new_direction = 'left' if _arguments['cell_id'] == 'left_cell' else 'right'
    else:
        _new_direction = _arguments['direction']
    _time_point_roi = roi(
        _length_x=_arguments['length_x'] * CELL_DIAMETER,
        _length_y=_arguments['length_y'] * CELL_DIAMETER,
        _offset_x=_arguments['offset_x'] * CELL_DIAMETER,
        _offset_y=_arguments['offset_y'] * CELL_DIAMETER,
        _cell_coordinates=_time_point_properties[_arguments['cell_id']]['coordinates'],
        _cell_diameter=_time_point_properties[_arguments['cell_id']]['diameter'],
        _direction=_new_direction
    )
    if _time_point_roi in _time_point_fibers_densities:
        return _arguments, _time_point_fibers_densities[_time_point_roi]
    else:
        if 'print' in _arguments and _arguments['print']:
            print('Computing:', _arguments['simulation'], _arguments['cell_id'], 'roi', _time_point_roi, 'direction', _arguments['direction'], 'tp',
                  _arguments['time_point'])
        _roi_fibers_density = roi_fibers_density(_arguments['simulation'], _arguments['time_point'], _time_point_roi)
        if 'save' in _arguments and _arguments['save']:
            _time_point_fibers_densities[_time_point_roi] = _roi_fibers_density
            save.fibers_densities(_arguments['simulation'], _arguments['time_point'], _time_point_fibers_densities)
        return _arguments, _roi_fibers_density


def roi_fibers_density_by_time(_arguments):
    _simulation_properties = properties(_arguments['simulation'])
    _rois = []
    for _time_point in range(min(_arguments['time_points'], len(_simulation_properties['time_points']))):
        _arguments['time_point'] = _time_point
        _rois.append(roi_fibers_density_time_point(_arguments)[1])

    return _arguments, _rois


def cells_distance(_properties):
    _c1_x = _properties['time_points'][0]['left_cell']['coordinates']['x']
    _c2_x = _properties['time_points'][0]['right_cell']['coordinates']['x']

    return round(abs(_c2_x - _c1_x) / CELL_DIAMETER, 10)


def is_in_region(_x, _y, _roi):
    return _roi[0] - EPSILON <= _x <= _roi[2] + EPSILON and _roi[1] - EPSILON <= _y <= _roi[3] + EPSILON


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

    return None


def get_expanded_roi(_roi):
    _expand_factor = 0.1
    return [_roi[0] - _expand_factor, _roi[1] - _expand_factor, _roi[2] + _expand_factor, _roi[3] + _expand_factor]


def get_master_roi(_rois):
    _x1, _y1, _x2, _y2 = -NET_DIMENSIONS['x1'], -NET_DIMENSIONS['y1'], -NET_DIMENSIONS['x2'], -NET_DIMENSIONS['y2']
    for _roi in _rois:
        _x1, _y1, _x2, _y2 = min(_x1, _roi[0]), min(_y1, _roi[1]), max(_x2, _roi[2]), max(_y2, _roi[3])

    return get_expanded_roi([_x1, _y1, _x2, _y2])


def elements_and_intersections_by_rois(_elements, _intersections, _tp0_fiber_lengths, _tp_fiber_lengths, _roi):
    _filtered_elements, _filtered_tp0_fibers_lengths, _filtered_tp_fibers_lengths = [], [], []
    for _element_index, _element in enumerate(_elements):
        _i1_id, _i2_id = _element
        _x1, _y1 = _intersections[_i1_id]
        _x2, _y2 = _intersections[_i2_id]
        # both intersections in region
        if is_in_region(_x1, _y1, _roi) and is_in_region(_x2, _y2, _roi):
            _filtered_elements.append(_element)
            if _tp0_fiber_lengths is not None:
                _filtered_tp0_fibers_lengths.append(_tp0_fiber_lengths[_element_index])
            if _tp_fiber_lengths is not None:
                _filtered_tp_fibers_lengths.append(_tp_fiber_lengths[_element_index])

    return _filtered_elements, _filtered_tp0_fibers_lengths, _filtered_tp_fibers_lengths
