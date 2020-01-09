from libs.experiments import load, paths, config, compute


def main():
    for _experiment in config.PAIRS:
        for _series in paths.folders(paths.fibers_density(_experiment)):
            _series_id = str(_series.split()[1])
            for _group in paths.folders(paths.fibers_density(_experiment, _series)):
                _group_id = _group.split('cells_')[1]
                for _z_group in paths.folders(paths.fibers_density(_experiment, _series, _group)):
                    _z_group_id = _z_group.split('cells_')[1]
                    _time_points_amount = len(
                        paths.text_files(paths.fibers_density(_experiment, _series, _group, _z_group))
                    )
                    _cell_coordinates_tracked_file = 'series_' + _series_id + '.txt'
                    _cell_coordinates_tracked = load.cell_coordinates_tracked_series_file_data(
                        _experiment, _cell_coordinates_tracked_file
                    )
                    _cell_1_id, _cell_2_id = _group_id.split('_')
                    _cell_1_coordinates = _cell_coordinates_tracked[int(_cell_1_id) - 1]
                    _cell_2_coordinates = _cell_coordinates_tracked[int(_cell_2_id) - 1]
                    _distance = compute.cells_distance_in_cell_size(
                        _experiment, _series, _cell_1_coordinates, _cell_2_coordinates
                    )

                    print(_experiment, _series_id, _group_id, _z_group_id, _time_points_amount, _distance, sep='\t')


if __name__ == '__main__':
    main()
