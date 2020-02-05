from libs.experiments import load, filtering


def main():
    _experiments = load.experiments_groups_as_tuples(['SN16', 'SN41'])
    _experiments = filtering.by_real_cells(_experiments, _real_cells=False)
    for _tuple in _experiments:
        _experiment, _series_id, _group = _tuple
        _properties = load.group_properties(_experiment, _series_id, _group)
        print(_properties['band'])


if __name__ == '__main__':
    main()
