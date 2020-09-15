import csv
import math
import os

from multiprocess.pool import Pool
from tqdm import tqdm

from libs import compute_lib
from libs.config_lib import CPUS_TO_USE
from libs.simulations import load, filtering, config, compute, paths, organize

PAIR_DISTANCE = [4, 5, 7, 9]
CAUSALITY = False
DOMINANT_PASSIVE = False
TIME_POINTS = math.inf


def main():
    _simulations = load.structured()
    _simulations = filtering.by_categories(
        _simulations,
        _is_single_cell=False,
        _is_heterogeneity=None,
        _is_low_connectivity=False,
        _is_causality=CAUSALITY,
        _is_dominant_passive=DOMINANT_PASSIVE
    )
    _simulations = filtering.by_pair_distances(_simulations, _distances=PAIR_DISTANCE)
    print('Total simulations:', len(_simulations))

    _arguments = []
    for _simulation in _simulations:
        for _cell_id in ['left_cell', 'right_cell']:
            _arguments.append({
                'simulation': _simulation,
                'length_x': config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                'length_y': config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                'offset_x': 0,
                'offset_y': 0,
                'cell_id': _cell_id,
                'direction': 'inside',
                'time_points': TIME_POINTS
            })

    _fiber_densities = {}
    with Pool(CPUS_TO_USE) as _p:
        for _keys, _value in tqdm(
                _p.imap_unordered(compute.window_fiber_density_by_time, _arguments),
                total=len(_arguments), desc='Computing windows & fiber densities'):
            _fiber_densities[(_keys['simulation'], _keys['cell_id'])] = _value
        _p.close()
        _p.join()

    _simulations_by_heterogeneity = organize.by_heterogeneity(_simulations)

    _headers = [
        'time_point',
        'simulation',
        'pair_distance_in_cell_diameter',
        'heterogeneity_std',
        'left_cell_fiber_density',
        'left_cell_fiber_density_z_score',
        'right_cell_fiber_density',
        'right_cell_fiber_density_z_score',
    ]

    _csv_path = os.path.join(paths.OUTPUTS, 'simulations_density_cell_pairs.csv')
    with open(_csv_path, 'w', newline='') as _csv_file:
        _csv_writer = csv.writer(_csv_file)
        _csv_writer.writerow(_headers)
        for _simulation_std in _simulations_by_heterogeneity:
            print('STD:', _simulation_std)
            _std_simulations = _simulations_by_heterogeneity[_simulation_std]
            for _simulation in tqdm(_std_simulations, desc='Simulations loop'):
                _properties = load.properties(_simulation)
                _pair_distance = compute.pair_distance(_properties)
                _average, _std = load.normalization(_simulation).values()
                _left_cell_fiber_densities = _fiber_densities[(_simulation, 'left_cell')]
                _right_cell_fiber_densities = _fiber_densities[(_simulation, 'right_cell')]
                for _time_point, (_left_cell_fiber_density, _right_cell_fiber_density) in \
                        enumerate(zip(_left_cell_fiber_densities, _right_cell_fiber_densities)):
                    _left_cell_fiber_density_z_score = compute_lib.z_score(_left_cell_fiber_density, _average, _std)
                    _right_cell_fiber_density_z_score = compute_lib.z_score(_right_cell_fiber_density, _average, _std)
                    _csv_writer.writerow([
                        _time_point,
                        _simulation,
                        _pair_distance,
                        _simulation_std,
                        _left_cell_fiber_density,
                        _left_cell_fiber_density_z_score,
                        _right_cell_fiber_density,
                        _right_cell_fiber_density_z_score
                    ])


if __name__ == '__main__':
    main()
