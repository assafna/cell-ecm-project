from libs import save_lib
from libs.experiments import paths, load

if __name__ == '__main__':
    with open('bands.txt', 'r') as _f:
        _lines = _f.readlines()
    for _experiment in ['SN16']:
        for _tuple in load.experiment_groups_as_tuples(_experiment):
            _experiment, _series_id, _group = _tuple
            _group_proproties = load.group_properties(_experiment, _series_id, _group)
            for _line in _lines:
                _series, _cells, _value = _line.split('\t')
                if _series == str(_series_id) and _group == _cells:
                    _group_proproties['band'] = True if _value == 'TRUE\n' else False
                    print(_experiment, _series_id, _group, _group_proproties['band'])
                    _path = paths.group_properties(_experiment, 'Series ' + str(_series_id), _group)
                    save_lib.to_json(_group_proproties, _path)
