import os

from libs import save_lib
from libs.experiments import load, paths


def main():
    _experiment = 'DeadDead_201208'
    _band = False

    _tuples = load.experiment_groups_as_tuples(_experiment)
    print('Total tuples:', len(_tuples))

    for _tuple in _tuples:
        print(_tuple)
        _experiment, _series_id, _group = _tuple
        _properties = load.group_properties(_experiment, _series_id, _group)
        _properties['band'] = _band

        _group_structured_path = paths.structured(_experiment, _series_id, _group)
        _properties_json_path = os.path.join(_group_structured_path, 'properties.json')
        save_lib.to_json(_properties, _properties_json_path)


if __name__ == '__main__':
    main()
