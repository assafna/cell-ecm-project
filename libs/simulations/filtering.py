from libs.simulations import load
from libs.simulations.compute import cells_distance


def by_time_points_amount(_simulations, _time_points, _exactly=False):
    if _exactly:
        return [_simulation for _simulation in _simulations if
                len(load.properties(_simulation)['time_points']) == _time_points]
    else:
        return [_simulation for _simulation in _simulations if
                len(load.properties(_simulation)['time_points']) >= _time_points]


def by_distance(_simulations, _distance):
    return [_simulation for _simulation in _simulations if
            cells_distance(load.properties(_simulation)) == _distance]


def by_distances(_simulations, _distances):
    return [_simulation for _simulation in _simulations if
            cells_distance(load.properties(_simulation)) in _distances]


def by_heterogeneity(_simulations, _std):
    return [_simulation for _simulation in _simulations if
            (_std == 0 and 'heterogeneity' not in _simulation)
            or (_std == 0.25 and 'std_025' in _simulation)
            or (_std == 0.5 and 'std' not in _simulation)
            or (_std == 0.75 and 'std_075' in _simulation)]


def is_single_cell(_simulation):
    return 'single_cell' in _simulation


def is_heterogeneity(_simulation):
    return 'heterogeneity' in _simulation


def is_low_connectivity(_simulation):
    return 'low_connectivity' in _simulation


def is_causality(_simulation):
    return 'causality' in _simulation


def is_dominant_passive(_simulation):
    return 'dominant_passive' in _simulation


def by_categories(_simulations, _is_single_cell=None, _is_heterogeneity=None, _is_low_connectivity=None,
                  _is_causality=None, _is_dominant_passive=None):
    return [_simulation for _simulation in _simulations if
            (_is_single_cell is None or _is_single_cell == is_single_cell(_simulation)) and
            (_is_heterogeneity is None or _is_heterogeneity == is_heterogeneity(_simulation)) and
            (_is_low_connectivity is None or _is_low_connectivity == is_low_connectivity(_simulation)) and
            (_is_causality is None or _is_causality == is_causality(_simulation)) and
            (_is_dominant_passive is None or _is_dominant_passive == is_dominant_passive(_simulation))
            ]
