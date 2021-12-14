import csv
import os

from tqdm import tqdm

from libs import compute_lib
from libs.experiments import config, filtering, load, compute, paths

PAIR_DISTANCE_RANGE = [4, 10]
BAND = None
REAL = None
STATIC = False
REAL_FAKE = False
OFFSETS_Y = [0, 0.5]


def main():
    _experiments = config.all_experiments()
    _experiments = filtering.by_categories(
        _experiments=_experiments,
        _is_single_cell=False,
        _is_high_temporal_resolution=None,
        _is_bleb=False,
        _is_dead_dead=False,
        _is_live_dead=False,
        _is_bead=False,
        _is_metastasis=False
    )

    _tuples = load.experiments_groups_as_tuples(_experiments)
    _tuples = filtering.by_pair_distance_range(_tuples, PAIR_DISTANCE_RANGE)
    _tuples = filtering.by_real_pairs(_tuples, _real_pairs=REAL)
    _tuples = filtering.by_fake_static_pairs(_tuples, _fake_static_pairs=STATIC)
    _tuples = filtering.by_real_fake_pairs(_tuples, _real_fake_pairs=REAL_FAKE)
    print('Total tuples:', len(_tuples))

    _arguments = []
    for _tuple in tqdm(_tuples, desc='Setting windows to compute'):
        _experiment, _series_id, _group = _tuple
        _latest_time_frame = compute.latest_time_frame_before_overlapping(_experiment, _series_id, _group, _offset_x=0)
        for _cell_id in ['left_cell', 'right_cell']:
            for _offset_y in OFFSETS_Y:
                _arguments.append({
                    'experiment': _experiment,
                    'series_id': _series_id,
                    'group': _group,
                    'length_x': config.QUANTIFICATION_WINDOW_LENGTH_IN_CELL_DIAMETER,
                    'length_y': config.QUANTIFICATION_WINDOW_HEIGHT_IN_CELL_DIAMETER,
                    'length_z': config.QUANTIFICATION_WINDOW_WIDTH_IN_CELL_DIAMETER,
                    'offset_x': 0,
                    'offset_y': _offset_y,
                    'offset_z': 0,
                    'cell_id': _cell_id,
                    'direction': 'inside',
                    'time_points': _latest_time_frame
                })

    _windows_dictionary, _windows_to_compute = \
        compute.windows(_arguments, _keys=['experiment', 'series_id', 'group', 'cell_id', 'offset_y'])
    _fiber_densities = compute.fiber_densities(_windows_to_compute)

    _experiments_fiber_densities = {
        _key: [_fiber_densities[_tuple] for _tuple in _windows_dictionary[_key]]
        for _key in _windows_dictionary
    }

    _headers = [
        'time_frame',
        'experiment',
        'series_id',
        'group',
        'left_cell_id',
        'right_cell_id',
        'band',
        'fake_following',
        'fake_static',
        'pair_distance_in_cell_diameter',
        'offset_z',
        'left_cell_fiber_density',
        'left_cell_fiber_density_z_score',
        'left_cell_out_of_boundaries',
        'right_cell_fiber_density',
        'right_cell_fiber_density_z_score',
        'right_cell_out_of_boundaries'
    ]

    _csv_path = os.path.join(paths.OUTPUTS, 'experiments_density_cell_pairs.csv')
    with open(_csv_path, 'w', newline='') as _csv_file:
        _csv_writer = csv.writer(_csv_file)
        _csv_writer.writerow(_headers)
        for _tuple in tqdm(_tuples, desc='Main loop'):
            _experiment, _series_id, _group = _tuple
            _group_properties = load.group_properties(_experiment, _series_id, _group)
            _left_cell_id, _right_cell_id = _group_properties['cells_ids'].values()
            _band = _group_properties['band']
            _fake_following, _fake_static = _group_properties['fake'], _group_properties['static']
            _average, _std = load.normalization_series_file_data(_experiment, _series_id).values()
            for _offset_y in OFFSETS_Y:
                _left_cell_fiber_densities = \
                    _experiments_fiber_densities[(_experiment, _series_id, _group, 'left_cell', _offset_y)]
                _right_cell_fiber_densities = \
                    _experiments_fiber_densities[(_experiment, _series_id, _group, 'right_cell', _offset_y)]
                _left_cell_fiber_densities = compute.remove_blacklist(
                    _experiment,
                    _series_id,
                    _left_cell_id,
                    _left_cell_fiber_densities
                )
                _right_cell_fiber_densities = compute.remove_blacklist(
                    _experiment,
                    _series_id,
                    _right_cell_id,
                    _right_cell_fiber_densities
                )
                for _time_frame, (_left_cell_fiber_density, _right_cell_fiber_density) in \
                        enumerate(zip(_left_cell_fiber_densities, _right_cell_fiber_densities)):
                    _pair_distance = compute.pair_distance_in_cell_size_time_frame(
                        _experiment, _series_id, _group, _time_frame
                    )
                    _left_cell_fiber_density, _left_cell_out_of_boundaries = _left_cell_fiber_density
                    _right_cell_fiber_density, _right_cell_out_of_boundaries = _right_cell_fiber_density
                    _left_cell_fiber_density_z_score = compute_lib.z_score(_left_cell_fiber_density, _average, _std)
                    _right_cell_fiber_density_z_score = compute_lib.z_score(_right_cell_fiber_density, _average, _std)
                    _csv_writer.writerow([
                        _time_frame,
                        _experiment,
                        _series_id,
                        _group,
                        _left_cell_id,
                        _right_cell_id,
                        _band,
                        _fake_following,
                        _fake_static,
                        _pair_distance,
                        _offset_y,
                        _left_cell_fiber_density,
                        _left_cell_fiber_density_z_score,
                        _left_cell_out_of_boundaries,
                        _right_cell_fiber_density,
                        _right_cell_fiber_density_z_score,
                        _right_cell_out_of_boundaries
                    ])


if __name__ == '__main__':
    main()
