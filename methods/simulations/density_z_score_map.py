import os
from itertools import product
from multiprocessing.pool import Pool

import numpy as np
import plotly.graph_objs as go
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import compute, paths
from libs.simulations import config
from libs.simulations import filtering
from libs.simulations import load
from plotting import save

OFFSET_Y = 0
TIME_POINT = 50
PAIR_DISTANCE = 7


def compute_fiber_densities(_simulations, _offsets_x, _time_point):
    _arguments = []
    for _simulation in _simulations:
        for _offset_x in _offsets_x:
            _arguments.append({
                'simulation': _simulation,
                'length_x': config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': _offset_x,
                'offset_y': OFFSET_Y,
                'cell_id': 'left_cell',
                'direction': 'inside',
                'time_point': _time_point
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_time_point, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['offset_x'])] = _value
        _p.close()
        _p.join()

    return _fiber_densities


def main():
    _simulations = load.structured()
    _simulations = filtering.by_time_points_amount(_simulations, _time_points=TIME_POINT)
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=False,
        _is_low_connectivity=False,
        _is_causality=False,
        _is_dominant_passive=False,
        _is_fibrin=False
    )
    _simulations = filtering.by_pair_distance(_simulations, _distance=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    # begin
    print('begin')
    _time_point = 0
    _offsets_x = [0, 1, 2, 3, 4, 5, 5]
    _fiber_densities = compute_fiber_densities(_simulations, _offsets_x, _time_point)
    for _offset_x in _offsets_x:
        _normalized_fiber_densities = []
        for _simulation in _simulations:
            _normalization = load.normalization(_simulation)
            _fiber_density = _fiber_densities[(_simulation, _offset_x)]

            _normalized_fiber_density = compute_lib.z_score(
                _fiber_density,
                _normalization['average'],
                _normalization['std']
            )
            _normalized_fiber_densities.append(_normalized_fiber_density)
        print(_offset_x, np.mean(_normalized_fiber_densities))

    # middle
    print('middle')
    _time_point = 25
    _offsets_x = [0, 1, 2, 3, 4, 5, 5, 5.25]
    _fiber_densities = compute_fiber_densities(_simulations, _offsets_x, _time_point)
    for _offset_x in _offsets_x:
        _normalized_fiber_densities = []
        for _simulation in _simulations:
            _normalization = load.normalization(_simulation)
            _fiber_density = _fiber_densities[(_simulation, _offset_x)]

            _normalized_fiber_density = compute_lib.z_score(
                _fiber_density,
                _normalization['average'],
                _normalization['std']
            )
            _normalized_fiber_densities.append(_normalized_fiber_density)
        print(_offset_x, np.mean(_normalized_fiber_densities))

    # end
    print('end')
    _time_point = 50
    _offsets_x = [0, 1, 2, 3, 4, 5, 5, 5.5]
    _fiber_densities = compute_fiber_densities(_simulations, _offsets_x, _time_point)
    for _offset_x in _offsets_x:
        _normalized_fiber_densities = []
        for _simulation in _simulations:
            _normalization = load.normalization(_simulation)
            _fiber_density = _fiber_densities[(_simulation, _offset_x)]

            _normalized_fiber_density = compute_lib.z_score(
                _fiber_density,
                _normalization['average'],
                _normalization['std']
            )
            _normalized_fiber_densities.append(_normalized_fiber_density)
        print(_offset_x, np.mean(_normalized_fiber_densities))


if __name__ == '__main__':
    main()
