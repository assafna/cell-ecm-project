import load_simulation
from computations import cells_distance


def by_time_points_amount(_simulations, _time_points, _exactly=False):
    if _exactly:
        return [_simulation for _simulation in _simulations if
                len(load_simulation.properties(_simulation)['time_points']) == _time_points]
    else:
        return [_simulation for _simulation in _simulations if
                len(load_simulation.properties(_simulation)['time_points']) >= _time_points]


def by_distance(_simulations, _distance):
    return [_simulation for _simulation in _simulations if
            cells_distance(load_simulation.properties(_simulation)) == _distance]
