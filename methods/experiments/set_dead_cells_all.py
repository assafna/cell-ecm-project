from libs.experiments import load, save


def main():
    _experiment = 'DeadDead_201208'
    _tuples = load.experiment_groups_as_tuples(_experiment)

    _cell_ids = {}
    for _tuple in _tuples:
        _experiment, _series_id, _group = _tuple
        _group_split = _group.split('_')
        _cell_1_id = int(_group_split[1])
        _cell_2_id = int(_group_split[2])

        if _series_id not in _cell_ids:
            _cell_ids[_series_id] = [_cell_1_id, _cell_2_id]
        else:
            if _cell_1_id not in _cell_ids[_series_id]:
                _cell_ids[_series_id].append(_cell_1_id)
            if _cell_2_id not in _cell_ids[_series_id]:
                _cell_ids[_series_id].append(_cell_2_id)

    for _series_id in _cell_ids.keys():
        _properties = load.image_properties(_experiment, _series_id)
        _properties['dead_cell_ids'] = _cell_ids[_series_id]
        save.image_properties(_experiment, _series_id, _properties)
        print(_experiment, _series_id, sep='\t')


if __name__ == '__main__':
    main()
