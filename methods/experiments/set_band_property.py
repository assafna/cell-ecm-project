import os

from multiprocess.pool import Pool

from libs import save_lib
from libs.config_lib import CPUS_TO_USE
from libs.experiments import compute, load, paths, config
from libs.experiments.config import AVERAGE_CELL_DIAMETER_IN_MICRONS


def process_group(_experiment, _series_id, _group, _overwrite=False):
    _group_properties = load.group_properties(_experiment, _series_id, _group)
    if not _overwrite and _group_properties['band'] is not None:
        return

    if 'static' in _group:
        return
    elif 'fake' in _group:
        _group_real = 'cells_' + _group.split('fake_')[1]
    else:
        _group_real = _group

    _group_real_properties = load.group_properties(_experiment, _series_id, _group_real)
    _time_points_amount = len(_group_real_properties['time_points'])

    _fibers_densities = []
    for _time_point in [0, int(round(_time_points_amount / 2)), _time_points_amount - 1]:
        _left_cell_coordinates = _group_real_properties['time_points'][_time_point]['left_cell']['coordinates']
        _right_cell_coordinates = _group_real_properties['time_points'][_time_point]['right_cell']['coordinates']
        _cell_diameter_x = AVERAGE_CELL_DIAMETER_IN_MICRONS / _group_real_properties['time_points'][_time_point]['resolutions']['x']
        _cell_diameter_y = AVERAGE_CELL_DIAMETER_IN_MICRONS / _group_real_properties['time_points'][_time_point]['resolutions']['y']
        _cell_diameter_z = AVERAGE_CELL_DIAMETER_IN_MICRONS / _group_real_properties['time_points'][_time_point]['resolutions']['z']
        _x1 = (_left_cell_coordinates['x'] + _right_cell_coordinates['x']) / 2 - _cell_diameter_x / 2
        _x2 = _x1 + _cell_diameter_x
        _y1 = (_left_cell_coordinates['y'] + _right_cell_coordinates['y']) / 2 - _cell_diameter_y / 2
        _y2 = _y1 + _cell_diameter_y
        _z1 = (_left_cell_coordinates['z'] + _right_cell_coordinates['z']) / 2 - _cell_diameter_z / 2
        _z2 = _z1 + _cell_diameter_z
        _roi = [_x1, _y1, _z1, _x2, _y2, _z2]
        _fibers_density = compute.roi_fibers_density(
            _experiment=_experiment,
            _series=_series_id,
            _group=_group_real,
            _time_point=_time_point,
            _roi=[int(round(_value)) for _value in _roi]
        )
        _fibers_densities.append(_fibers_density)

    if _fibers_densities[0] < _fibers_densities[1] < _fibers_densities[2]:
        _band = True
    else:
        _band = False

    print(_experiment, _series_id, _group, _band, sep='\t')

    _group_properties['band'] = _band
    _group_structured_path = paths.structured(_experiment, 'Series ' + str(_series_id), _group)
    _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
    save_lib.to_json(_group_properties, _properties_json_path)


def process_experiment(_experiment, _overwrite=False):
    _arguments = []
    for _tuple in load.experiment_groups_as_tuples(_experiment):
        _experiment, _series_id, _group = _tuple
        _arguments.append((_experiment, _series_id, _group, _overwrite))

    _p = Pool(CPUS_TO_USE)
    _p.starmap(process_group, _arguments)
    _p.close()


def process_experiments(_experiments, _overwrite=False):
    for _experiment in _experiments:
        process_experiment(_experiment, _overwrite)


def process_all_experiments(_overwrite=False):
    process_experiments(config.PAIRS, _overwrite)


if __name__ == '__main__':
    process_all_experiments()
