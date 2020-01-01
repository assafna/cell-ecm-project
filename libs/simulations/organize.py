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
