from libs.experiments import load, paths, config, compute


def process_group(_experiment, _series_id, _group):
    _properties = load.group_properties(_experiment, _series_id, _group)
    _cell_1_coordinates = _properties['time_points'][0]['left_cell']['coordinates'].values()
    _cell_2_coordinates = _properties['time_points'][0]['right_cell']['coordinates'].values()
    _distance = compute.cells_distance_in_cell_size(
        _experiment, _series_id, [_cell_1_coordinates], [_cell_2_coordinates]
    )

    print(_experiment, _series_id, _group, _distance, _properties['band'], _properties['fake'], _properties['static'], len(_properties['time_points']), sep='\t')


def process_experiment(_experiment):
    for _tuple in load.experiment_groups_as_tuples(_experiment):
        _experiment, _series_id, _group = _tuple
        process_group(_experiment, _series_id, _group)


def process_experiments(_experiments):
    for _experiment in _experiments:
        process_experiment(_experiment)


def process_all_experiments():
    # TODO: handle single cell
    process_experiments(config.PAIRS)


if __name__ == '__main__':
    process_all_experiments()
