from libs.simulations import compute, load


def by_distances(_simulations, _reverse=False):
    _organized_simulations = {}
    for _simulation in _simulations:
        _simulation_properties = load.properties(_simulation)
        _cells_distance = compute.cells_distance(_simulation_properties)
        if _cells_distance in _organized_simulations:
            _organized_simulations[_cells_distance].append(_simulation)
        else:
            _organized_simulations[_cells_distance] = [_simulation]

    return {_distance: _organized_simulations[_distance] for
            _distance in sorted(_organized_simulations.keys(), reverse=_reverse)}


def by_heterogeneity(_simulations):
    _organized_simulations = {}
    for _simulation in _simulations:
        if '_heterogeneity' not in _simulation:
            _std = 0
        elif '_std_' not in _simulation:
            _std = 0.5
        elif '_std_025' in _simulation:
            _std = 0.25
        elif '_std_075' in _simulation:
            _std = 0.75
        else:
            raise Exception('No such STD')
        if _std in _organized_simulations:
            _organized_simulations[_std].append(_simulation)
        else:
            _organized_simulations[_std] = [_simulation]

    return _organized_simulations
