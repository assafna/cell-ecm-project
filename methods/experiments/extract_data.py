from libs.experiments import load, paths, config, compute


def main():
    for _tuple in load.experiment_groups_as_tuples('SN41'):
        _experiment, _series_id, _group = _tuple
        _cell_coordinates_tracked_file = 'series_' + str(_series_id) + '.txt'
        _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(
            _experiment, _cell_coordinates_tracked_file
        )
        _, _cell_1_id, _cell_2_id = _group.split('_')
        _cell_1_coordinates = _cell_coordinates_tracked[int(_cell_1_id)]
        _cell_2_coordinates = _cell_coordinates_tracked[int(_cell_2_id)]
        _distance = compute.cells_distance_in_cell_size(
            _experiment, _series_id, _cell_1_coordinates, _cell_2_coordinates
        )
        _properties = load.group_properties(_experiment, _series_id, _group)

        print(_experiment, _series_id, _group, _distance, _properties['band'], len(_properties['time_points']), sep='\t')


if __name__ == '__main__':
    main()
