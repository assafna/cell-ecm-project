from libs.simulations import load
from libs.simulations.compute import pair_distance


def by_time_points_amount(_simulations, _time_points, _exactly=False):
    if _exactly:
        return [_simulation for _simulation in _simulations if
                len(load.properties(_simulation)['time_points']) == _time_points]
    else:
        return [_simulation for _simulation in _simulations if
                len(load.properties(_simulation)['time_points']) >= _time_points]


def by_pair_distance(_simulations, _distance):
    return [_simulation for _simulation in _simulations if
            pair_distance(load.properties(_simulation)) == _distance]


def by_pair_distances(_simulations, _distances):
    return [_simulation for _simulation in _simulations if
            pair_distance(load.properties(_simulation)) in _distances]


def by_heterogeneity(_simulations, _std):
    return [_simulation for _simulation in _simulations if
            (_std == 0 and 'heterogeneity' not in _simulation)
            or (_std == 0.25 and 'std_025' in _simulation)
            or (_std == 0.5 and 'heterogeneity' in _simulation and 'std' not in _simulation)
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


def by_causality(_simulations, _alpha=None, _beta=None):
    if _alpha is None and _beta is None:
        return _simulations

    _alpha_str = ''
    if _alpha is not None:
        if _alpha == 0:
            _alpha_str = 'alpha_000'
        elif _alpha == 0.25:
            _alpha_str = 'alpha_025'
        elif _alpha == 0.5:
            _alpha_str = 'alpha_050'
        elif _alpha == 0.75:
            _alpha_str = 'alpha_075'
        elif _alpha == 1:
            _alpha_str = 'alpha_100'

    _beta_str = ''
    if _beta is not None:
        if _beta == 1:
            _beta_str = 'beta_100'
        elif _beta == 1.05:
            _beta_str = 'beta_105'
        elif _beta == 1.1:
            _beta_str = 'beta_110'
        elif _beta == 1.2:
            _beta_str = 'beta_120'

    _simulations_filtered = []
    for _simulation in _simulations:
        if _alpha is not None:
            if _alpha > 0:
                if _beta is not None:
                    if _alpha_str in _simulation and _beta_str in _simulation:
                        _simulations_filtered.append(_simulation)
                else:
                    if _alpha_str in _simulation:
                        _simulations_filtered.append(_simulation)
            # alpha is 0, beta does not matter
            else:
                if _alpha_str in _simulation:
                    _simulations_filtered.append(_simulation)
        # alpha does not matter
        else:
            if _beta_str in _simulation:
                _simulations_filtered.append(_simulation)

    return _simulations_filtered
